//! Format Anagène « .edi » — portage fidèle de `geniegen/core/edi.py`.
//!
//! Le `.edi` est un format TEXTE en **Latin-1** (Windows-1252), structuré par
//! des lignes de commentaire commençant par `;`. Il n'est pas propriétaire :
//! ce module est la couche d'interopérabilité qui permet à DZ-Gene de lire les
//! banques de séquences déjà utilisées dans les lycées.
//!
//! ```text
//! ; Anagène - Fenêtre Edition     <- en-tête du fichier
//! ; HBB-normal                    <- 1re ligne de méta = nom
//! ; Type 1                        <- 1=ADN, 2=protéine 3-lettres, 3=prot. 1-lettre
//! ; Dec  0                        <- décalage visuel (positions vides à gauche)
//! ; Allèle sain de la bêta-globine
//! ATGGTGCACCTGACTCCTGAGGAG
//! ;-                              <- séparateur de blocs
//! ```
//!
//! **Contrat de non-régression :** ce module doit lire *sans perte* les
//! fichiers de `banque-dz/`. Un test d'aller-retour (parse → write → parse)
//! le garantit.

use serde::{Deserialize, Serialize};
use std::fmt;

/// En-tête attendu en tête de fichier.
pub const HEADER: &str = "; Anagène - Fenêtre Edition";
/// Séparateur de blocs.
pub const SEPARATOR: &str = ";-";
/// Longueur de repliement des séquences à l'écriture.
pub const LINE_WRAP: usize = 68;
const EOL: &str = "\r\n";

/// Type de séquence, tel que codé dans le champ `Type` du `.edi`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SeqType {
    /// Acide nucléique (A, C, G, T ; ARN toléré avec U).
    Adn,
    /// Protéine, code à 3 lettres (MetAlaGlu…).
    Prot3,
    /// Protéine, code à 1 lettre (MAE…).
    Prot1,
}

impl SeqType {
    /// Vrai si le type code une protéine (3 ou 1 lettre).
    pub fn is_protein(self) -> bool {
        matches!(self, SeqType::Prot3 | SeqType::Prot1)
    }

    /// Code brut du champ `Type`.
    pub fn code(self) -> u8 {
        match self {
            SeqType::Adn => 1,
            SeqType::Prot3 => 2,
            SeqType::Prot1 => 3,
        }
    }

    /// Type depuis le code brut. Tout code inconnu retombe sur `Adn`
    /// (tolérance délibérée : mieux vaut afficher que refuser en classe).
    pub fn from_code(code: u8) -> Self {
        match code {
            2 => SeqType::Prot3,
            3 => SeqType::Prot1,
            _ => SeqType::Adn,
        }
    }

    /// Libellé lisible du type de séquence.
    pub fn label(self) -> &'static str {
        match self {
            SeqType::Adn => "ADN/ARN",
            SeqType::Prot3 => "Protéine (3 lettres)",
            SeqType::Prot1 => "Protéine (1 lettre)",
        }
    }
}

/// Table code 3 lettres → 1 lettre.
const AA3_TO_1: [(&str, char); 27] = [
    ("Ala", 'A'), ("Arg", 'R'), ("Asn", 'N'), ("Asp", 'D'), ("Cys", 'C'),
    ("Gln", 'Q'), ("Glu", 'E'), ("Gly", 'G'), ("His", 'H'), ("Ile", 'I'),
    ("Leu", 'L'), ("Lys", 'K'), ("Met", 'M'), ("Phe", 'F'), ("Pro", 'P'),
    ("Ser", 'S'), ("Thr", 'T'), ("Trp", 'W'), ("Tyr", 'Y'), ("Val", 'V'),
    ("Sec", 'U'), ("Pyl", 'O'), ("Asx", 'B'), ("Glx", 'Z'), ("Xaa", 'X'),
    ("Ter", '*'), ("End", '*'),
];

fn aa3_to_1(triplet: &str) -> Option<char> {
    AA3_TO_1
        .iter()
        .find(|(k, _)| k.eq_ignore_ascii_case(triplet))
        .map(|(_, v)| *v)
}

fn aa1_to_3(aa: char) -> &'static str {
    if aa == '*' {
        return "Ter";
    }
    AA3_TO_1
        .iter()
        .find(|(_, v)| *v == aa.to_ascii_uppercase())
        .map(|(k, _)| *k)
        .unwrap_or("Xaa")
}

/// Une séquence de la banque.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sequence {
    /// Nom affiché / interne.
    pub name: String,
    /// Type de la séquence (ADN, protéine 3 ou 1 lettre).
    pub seqtype: SeqType,
    /// Séquence **normalisée** : les protéines sont toujours ramenées au code
    /// 1 lettre, quel que soit le `Type` d'origine.
    pub residues: String,
    /// Décalage : positions vides insérées à gauche pour l'alignement visuel.
    pub dec: i32,
    /// Lignes de commentaire additionnelles.
    pub description: String,
}

impl Sequence {
    /// Construit une séquence minimale (dec=0, sans description).
    pub fn new(name: impl Into<String>, seqtype: SeqType, residues: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            seqtype,
            residues: residues.into(),
            dec: 0,
            description: String::new(),
        }
    }

    /// Nombre de résidus de la séquence.
    pub fn length(&self) -> usize {
        self.residues.len()
    }

    /// Rend une protéine au format 3 lettres (concaténé, sans espaces).
    pub fn as_three_letter(&self) -> Result<String, EdiError> {
        if !self.seqtype.is_protein() {
            return Err(EdiError::NotAProtein);
        }
        Ok(self.residues.chars().map(aa1_to_3).collect())
    }
}

/// Un fichier `.edi` = une liste ordonnée de séquences.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EdiFile {
    /// Séquences du fichier, dans l'ordre de lecture.
    pub sequences: Vec<Sequence>,
}

impl EdiFile {
    /// Nombre de séquences dans le fichier.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }
    /// Vrai si le fichier ne contient aucune séquence.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
    /// Itère sur les séquences du fichier.
    pub fn iter(&self) -> std::slice::Iter<'_, Sequence> {
        self.sequences.iter()
    }
}

/// Erreurs du module `.edi`.
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum EdiError {
    #[error("le format 3-lettres ne s'applique qu'aux protéines")]
    /// Opération réservée aux protéines, appliquée à un acide nucléique.
    NotAProtein,
}

// --------------------------------------------------------------------------
// Encodage Latin-1
// --------------------------------------------------------------------------

/// Décode des octets Latin-1 en `String`.
///
/// Latin-1 est une bijection octet → point de code U+0000..U+00FF : le décodage
/// ne peut **jamais** échouer. C'est précisément pourquoi Anagène l'utilise, et
/// pourquoi on ne doit surtout pas tenter un décodage UTF-8 ici (les « é » des
/// descriptions françaises feraient échouer la lecture de toute la banque).
pub fn decode_latin1(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| b as char).collect()
}

/// Encode une `String` en Latin-1. Les caractères hors Latin-1 (arabe…) sont
/// remplacés par `?` — comme le ferait Anagène lui-même.
pub fn encode_latin1(s: &str) -> Vec<u8> {
    s.chars()
        .map(|c| if (c as u32) < 256 { c as u8 } else { b'?' })
        .collect()
}

// --------------------------------------------------------------------------
// Normalisation des résidus
// --------------------------------------------------------------------------

/// Ne garde que lettres, `-` (gap) et `*` (stop) ; met en majuscules.
fn clean_nucleic(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphabetic() || *c == '-' || *c == '*')
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

/// Convertit une chaîne de triplets (`MetAlaGlu…`) en code 1 lettre.
///
/// Un triplet inconnu devient `X` — on avance quand même de 3 caractères pour
/// rester aligné : en classe, une banque légèrement corrompue doit s'ouvrir.
fn prot3_to_1(raw: &str) -> String {
    let compact: Vec<char> = raw
        .chars()
        .filter(|c| c.is_ascii_alphabetic() || *c == '*')
        .collect();

    compact
        .chunks(3)
        .map(|chunk| {
            let trip: String = chunk.iter().collect();
            aa3_to_1(&trip).unwrap_or('X')
        })
        .collect()
}

// --------------------------------------------------------------------------
// Lecture
// --------------------------------------------------------------------------

/// Extrait la valeur d'une méta-ligne `; Clé  valeur`.
fn meta_value<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let body = line.strip_prefix(';')?.trim();
    let rest = body.strip_prefix(key).or_else(|| {
        // comparaison insensible à la casse, sans allouer
        body.get(..key.len())
            .filter(|p| p.eq_ignore_ascii_case(key))
            .and_then(|_| body.get(key.len()..))
    })?;
    Some(rest.trim())
}

/// Parse le contenu textuel d'un fichier `.edi`.
///
/// Robuste par conception : une ligne mal formée est ignorée plutôt que de
/// faire échouer toute la banque.
#[allow(unused_assignments)] // les resets en fin de macro flush! sont volontaires
pub fn parse(text: &str) -> EdiFile {
    let mut sequences = Vec::new();

    // Bloc en cours de construction.
    let mut name: Option<String> = None;
    let mut raw_type: u8 = 1;
    let mut dec: i32 = 0;
    let mut description: Vec<String> = Vec::new();
    let mut residues = String::new();
    let mut in_block = false;

    // Ferme le bloc courant et l'ajoute au résultat.
    macro_rules! flush {
        () => {
            if in_block {
                if let Some(n) = name.take() {
                    let seqtype = SeqType::from_code(raw_type);
                    let normalized = if seqtype == SeqType::Prot3 {
                        prot3_to_1(&residues)
                    } else {
                        clean_nucleic(&residues)
                    };
                    sequences.push(Sequence {
                        name: n,
                        seqtype,
                        residues: normalized,
                        dec,
                        description: description.join("\n"),
                    });
                }
            }
            name = None;
            raw_type = 1;
            dec = 0;
            description.clear();
            residues.clear();
            in_block = false;
        };
    }

    for line in text.lines() {
        let trimmed = line.trim_end_matches(['\r', '\n']);

        if trimmed.trim() == SEPARATOR {
            flush!();
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix(';') {
            let body = rest.trim();

            // En-tête de fichier : ignoré.
            if body.starts_with("Anagène") || body.starts_with("Anag") {
                continue;
            }

            if let Some(v) = meta_value(trimmed, "Type") {
                raw_type = v.parse().unwrap_or(1);
                continue;
            }
            if let Some(v) = meta_value(trimmed, "Dec") {
                dec = v.parse().unwrap_or(0);
                continue;
            }

            // 1re méta non reconnue = nom ; les suivantes = description.
            if name.is_none() {
                name = Some(body.to_owned());
                in_block = true;
            } else if !body.is_empty() {
                description.push(body.to_owned());
            }
            continue;
        }

        // Ligne de séquence.
        if !trimmed.trim().is_empty() {
            residues.push_str(trimmed.trim());
            in_block = true;
        }
    }

    flush!();
    EdiFile { sequences }
}

/// Parse des octets Latin-1 bruts (le cas réel : lecture d'un fichier).
pub fn parse_bytes(bytes: &[u8]) -> EdiFile {
    parse(&decode_latin1(bytes))
}

// --------------------------------------------------------------------------
// Écriture
// --------------------------------------------------------------------------

impl fmt::Display for EdiFile {
    /// Sérialise au format `.edi` (CRLF, repliement à 68 caractères).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{HEADER}{EOL}")?;

        for seq in &self.sequences {
            write!(f, "; {}{EOL}", seq.name)?;
            write!(f, "; Type {}{EOL}", seq.seqtype.code())?;
            write!(f, "; Dec  {}{EOL}", seq.dec)?;

            for line in seq.description.lines().filter(|l| !l.trim().is_empty()) {
                write!(f, "; {line}{EOL}")?;
            }

            // Les protéines Type 2 sont ré-exportées en 3 lettres, fidèlement.
            let body = if seq.seqtype == SeqType::Prot3 {
                seq.as_three_letter().unwrap_or_else(|_| seq.residues.clone())
            } else {
                seq.residues.clone()
            };

            let chars: Vec<char> = body.chars().collect();
            for chunk in chars.chunks(LINE_WRAP) {
                let line: String = chunk.iter().collect();
                write!(f, "{line}{EOL}")?;
            }

            write!(f, "{SEPARATOR}{EOL}")?;
        }
        Ok(())
    }
}

/// Sérialise en octets Latin-1, prêts à écrire sur disque.
pub fn to_bytes(file: &EdiFile) -> Vec<u8> {
    encode_latin1(&file.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = concat!(
        "; Anagène - Fenêtre Edition\r\n",
        "; HBB-normal\r\n",
        "; Type 1\r\n",
        "; Dec  0\r\n",
        "; Allèle sain de la bêta-globine\r\n",
        "ATGGTGCACCTGACTCCTGAGGAG\r\n",
        ";-\r\n",
        "; HBB-drepano\r\n",
        "; Type 1\r\n",
        "; Dec  0\r\n",
        "ATGGTGCACCTGACTCCTGTGGAG\r\n",
        ";-\r\n",
    );

    #[test]
    fn parses_two_sequences() {
        let f = parse(SAMPLE);
        assert_eq!(f.len(), 2);
        assert_eq!(f.sequences[0].name, "HBB-normal");
        assert_eq!(f.sequences[0].seqtype, SeqType::Adn);
        assert_eq!(f.sequences[0].residues, "ATGGTGCACCTGACTCCTGAGGAG");
        assert_eq!(f.sequences[1].name, "HBB-drepano");
    }

    #[test]
    fn keeps_accented_description() {
        let f = parse(SAMPLE);
        assert!(f.sequences[0].description.contains("bêta-globine"));
    }

    #[test]
    fn latin1_never_fails_to_decode() {
        // 0xEA = « ê » en Latin-1 ; ce serait un octet UTF-8 invalide.
        let bytes = [b';', b' ', b'b', 0xEA, b't', b'a'];
        assert_eq!(decode_latin1(&bytes), "; bêta");
    }

    #[test]
    fn latin1_round_trip() {
        let s = "Allèle sain de la bêta-globine";
        assert_eq!(decode_latin1(&encode_latin1(s)), s);
    }

    #[test]
    fn three_letter_protein_is_normalised_to_one_letter() {
        let text = "; Anagène - Fenêtre Edition\r\n\
                    ; Prot-test\r\n; Type 2\r\n; Dec  0\r\nMetAlaGluTer\r\n;-\r\n";
        let f = parse(text);
        assert_eq!(f.sequences[0].residues, "MAE*");
        assert_eq!(f.sequences[0].seqtype, SeqType::Prot3);
    }

    #[test]
    fn three_letter_export_is_faithful() {
        let text = "; Anagène - Fenêtre Edition\r\n\
                    ; P\r\n; Type 2\r\n; Dec  0\r\nMetAlaGlu\r\n;-\r\n";
        let f = parse(text);
        assert_eq!(f.sequences[0].as_three_letter().unwrap(), "MetAlaGlu");
    }

    #[test]
    fn round_trip_preserves_content() {
        let original = parse(SAMPLE);
        let reparsed = parse(&original.to_string());
        assert_eq!(original, reparsed, "parse -> write -> parse doit être stable");
    }

    #[test]
    fn dec_offset_is_read() {
        let text = "; Anagène - Fenêtre Edition\r\n\
                    ; S\r\n; Type 1\r\n; Dec  12\r\nATGC\r\n;-\r\n";
        assert_eq!(parse(text).sequences[0].dec, 12);
    }

    #[test]
    fn long_sequence_is_wrapped_at_68() {
        let long = "A".repeat(150);
        let f = EdiFile { sequences: vec![Sequence::new("long", SeqType::Adn, &long)] };
        let out = f.to_string();
        let body: Vec<&str> = out
            .lines()
            .filter(|l| l.starts_with('A'))
            .collect();
        assert_eq!(body.len(), 3); // 68 + 68 + 14
        assert_eq!(body[0].len(), LINE_WRAP);
        // et la relecture doit rendre la séquence intacte
        assert_eq!(parse(&out).sequences[0].residues, long);
    }

    #[test]
    fn malformed_lines_do_not_crash() {
        let text = "; Anagène - Fenêtre Edition\r\n\
                    ; S\r\n; Type abc\r\n; Dec  xyz\r\nATGC\r\n;-\r\n";
        let f = parse(text);
        assert_eq!(f.len(), 1);
        assert_eq!(f.sequences[0].seqtype, SeqType::Adn); // repli sûr
        assert_eq!(f.sequences[0].dec, 0);
    }

    #[test]
    fn empty_input_yields_empty_file() {
        assert!(parse("").is_empty());
    }
}
