//! Moteur génétique — portage fidèle de `geniegen/core/genetics.py`.
//!
//! Opérations pédagogiques du menu « Traitement » d'Anagène :
//! transcription, traduction, ORF, comparaison simple.
//!
//! Rust pur : aucune allocation cachée, aucune dépendance scientifique externe.
//! Ces opérations doivent rester utilisables **hors-ligne, en salle de classe**.

use serde::{Deserialize, Serialize};

/// Code génétique standard (NCBI table 1).
///
/// Indexation : `base(pos1) * 16 + base(pos2) * 4 + base(pos3)`
/// avec l'ordre canonique NCBI `T/U = 0, C = 1, A = 2, G = 3`.
///
/// Cette table de 64 octets remplace le `dict` Python : lookup en O(1)
/// sans hachage ni allocation.
const AA_TABLE: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Codon d'initiation.
pub const START_CODON: &str = "AUG";
/// Codons stop.
pub const STOP_CODONS: [&str; 3] = ["UAA", "UAG", "UGA"];

/// Indice canonique d'une base (`T` et `U` sont équivalents).
#[inline]
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Traduit un codon en acide aminé (code 1 lettre).
///
/// Retourne `None` si le codon n'a pas exactement 3 bases valides
/// (ambiguïtés `N`, `R`, `Y`… incluses) — l'appelant décide alors du `X`.
///
/// ```
/// use dzgene_core::genetics::codon_to_aa;
/// assert_eq!(codon_to_aa("AUG"), Some('M'));
/// assert_eq!(codon_to_aa("ATG"), Some('M')); // T et U équivalents
/// assert_eq!(codon_to_aa("UAA"), Some('*'));
/// assert_eq!(codon_to_aa("NNN"), None);
/// ```
pub fn codon_to_aa(codon: &str) -> Option<char> {
    let b = codon.as_bytes();
    if b.len() != 3 {
        return None;
    }
    let i = base_index(b[0])? * 16 + base_index(b[1])? * 4 + base_index(b[2])?;
    Some(AA_TABLE[i] as char)
}

/// Complément d'une base d'ADN.
#[inline]
fn complement_dna_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other, // gaps `-`, ambiguïtés : laissés intacts
    }
}

/// ADN → ARN (remplace T par U).
pub fn to_rna(dna: &str) -> String {
    dna.to_ascii_uppercase().replace('T', "U")
}

/// ARN → ADN (remplace U par T).
pub fn to_dna(rna: &str) -> String {
    rna.to_ascii_uppercase().replace('U', "T")
}

/// Brin complémentaire (même sens de lecture).
pub fn complement(dna: &str) -> String {
    dna.bytes().map(|b| complement_dna_base(b) as char).collect()
}

/// Brin complémentaire inverse (lu 5'→3').
pub fn reverse_complement(dna: &str) -> String {
    dna.bytes()
        .rev()
        .map(|b| complement_dna_base(b) as char)
        .collect()
}

/// Transcription à partir du **brin matrice**.
///
/// Le transcrit est le complémentaire inverse du brin matrice : il a donc la
/// même séquence que le brin codant, avec U à la place de T.
pub fn transcribe(template_dna: &str) -> String {
    to_rna(&reverse_complement(template_dna))
}

/// Résultat d'une traduction.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Translation {
    /// Séquence protéique en code 1 lettre (`*` = STOP, `X` = codon illisible).
    pub protein: String,
    /// Codons ARN effectivement lus, dans l'ordre.
    pub codons: Vec<String>,
    /// Cadre de lecture utilisé (0, 1 ou 2).
    pub frame: usize,
    /// `true` si un codon STOP a été rencontré.
    pub stopped: bool,
}

/// Traduit une séquence nucléique en protéine.
///
/// * `frame`   — cadre de lecture 0/1/2.
/// * `to_stop` — arrête la lecture au premier STOP (le `*` reste inclus).
///
/// L'ARN est détecté automatiquement (présence de `U` sans `T`), comme en
/// Python ; `to_rna` est idempotent, donc la conversion est sûre dans tous les
/// cas.
pub fn translate(seq: &str, frame: usize, to_stop: bool) -> Translation {
    let rna = to_rna(seq);
    let bytes = rna.as_bytes();

    let mut protein = String::new();
    let mut codons = Vec::new();
    let mut stopped = false;

    let mut i = frame;
    while i + 3 <= bytes.len() {
        let codon = &rna[i..i + 3];
        codons.push(codon.to_owned());

        let aa = codon_to_aa(codon).unwrap_or('X');
        protein.push(aa);

        if aa == '*' {
            stopped = true;
            if to_stop {
                break;
            }
        }
        i += 3;
    }

    Translation { protein, codons, frame, stopped }
}

/// Un cadre ouvert de lecture.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Orf {
    /// Position nucléique de début (0-based, sur le `AUG`).
    pub start: usize,
    /// Position nucléique de fin (exclusive, après le codon STOP).
    pub end: usize,
    /// Protéine produite, **sans** le `*` final.
    pub protein: String,
}

/// Cherche le **premier** cadre ouvert de lecture (AUG → STOP).
///
/// Balaie les 3 cadres, de gauche à droite. Retourne `None` si aucun AUG
/// n'est suivi d'un STOP dans le même cadre.
pub fn find_orf(seq: &str) -> Option<Orf> {
    let rna = to_rna(seq);

    for frame in 0..3 {
        let mut i = frame;
        while i + 3 <= rna.len() {
            if &rna[i..i + 3] == START_CODON {
                let tr = translate(&rna[i..], 0, true);
                if tr.stopped {
                    let end = i + tr.protein.len() * 3;
                    return Some(Orf {
                        start: i,
                        end,
                        protein: tr.protein.trim_end_matches('*').to_owned(),
                    });
                }
            }
            i += 3;
        }
    }
    None
}

/// Résultat d'une « comparaison simple » (sans indel), à la manière d'Anagène.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SimpleComparison {
    /// Longueur comparée = min des deux longueurs.
    pub length: usize,
    /// Nombre de positions identiques.
    pub identical: usize,
    /// Positions (0-based) qui diffèrent.
    pub differences: Vec<usize>,
    /// Pourcentage d'identité sur la longueur comparée.
    pub identity_pct: f64,
}

impl SimpleComparison {
    /// Nombre de positions qui diffèrent entre les deux séquences.
    pub fn n_differences(&self) -> usize {
        self.differences.len()
    }
}

/// Compare deux séquences position par position, **sans décalage**.
pub fn compare_simple(a: &str, b: &str) -> SimpleComparison {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();

    let differences: Vec<usize> = a
        .bytes()
        .zip(b.bytes())
        .enumerate()
        .filter_map(|(i, (x, y))| (x != y).then_some(i))
        .collect();

    let length = a.len().min(b.len());
    let identical = length - differences.len();
    let identity_pct = if length == 0 {
        0.0
    } else {
        identical as f64 / length as f64 * 100.0
    };

    SimpleComparison { length, identical, differences, identity_pct }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_table_is_ncbi_table_1() {
        // Ancrages sur les codons que les élèves de 3AS doivent connaître.
        assert_eq!(codon_to_aa("AUG"), Some('M')); // initiation
        assert_eq!(codon_to_aa("UAA"), Some('*'));
        assert_eq!(codon_to_aa("UAG"), Some('*'));
        assert_eq!(codon_to_aa("UGA"), Some('*'));
        assert_eq!(codon_to_aa("UGG"), Some('W')); // seul codon du Trp
        assert_eq!(codon_to_aa("GGG"), Some('G'));
        assert_eq!(codon_to_aa("UUU"), Some('F'));
        // Dégénérescence : les 4 codons GC* codent tous l'alanine.
        for c in ["GCU", "GCC", "GCA", "GCG"] {
            assert_eq!(codon_to_aa(c), Some('A'), "codon {c}");
        }
    }

    #[test]
    fn t_and_u_are_interchangeable() {
        assert_eq!(codon_to_aa("ATG"), codon_to_aa("AUG"));
    }

    #[test]
    fn ambiguous_codons_have_no_aa() {
        assert_eq!(codon_to_aa("NNN"), None);
        assert_eq!(codon_to_aa("AU"), None);
        assert_eq!(codon_to_aa("AUGC"), None);
    }

    #[test]
    fn reverse_complement_is_an_involution() {
        let s = "ATGCGTTAAC";
        assert_eq!(reverse_complement(&reverse_complement(s)), s);
    }

    #[test]
    fn transcription_from_template_strand() {
        // Le brin matrice est lu 3'->5' par l'ARN polymerase, qui synthetise
        // l'ARNm par complementarite et en sens inverse. Donc :
        //   matrice CAT -> complement inverse ATG -> ARNm AUG (codon START).
        assert_eq!(transcribe("CAT"), "AUG");

        // Piege classique : la matrice TAC ne donne PAS AUG, mais GUA (Val).
        assert_eq!(transcribe("TAC"), "GUA");
    }

    #[test]
    fn translation_stops_at_stop_codon() {
        // AUG GCU UAA GGG : Met-Ala-STOP, puis Gly si on ne s'arrête pas.
        let full = translate("AUGGCUUAAGGG", 0, false);
        assert_eq!(full.protein, "MA*G");
        assert!(full.stopped);

        let short = translate("AUGGCUUAAGGG", 0, true);
        assert_eq!(short.protein, "MA*");
        assert_eq!(short.codons.len(), 3);
    }

    #[test]
    fn reading_frame_shifts_the_protein() {
        let f0 = translate("AAUGGCUUAA", 0, false);
        let f1 = translate("AAUGGCUUAA", 1, false);
        assert_ne!(f0.protein, f1.protein);
        assert_eq!(f1.protein, "MA*"); // le cadre +1 révèle l'ORF
    }

    #[test]
    fn trailing_incomplete_codon_is_ignored() {
        // 8 bases = 2 codons complets + 2 bases orphelines.
        let tr = translate("AUGGCUUA", 0, false);
        assert_eq!(tr.codons.len(), 2);
        assert_eq!(tr.protein, "MA");
    }

    #[test]
    fn find_orf_locates_first_aug_to_stop() {
        let orf = find_orf("CCAUGGCUUAAGG").expect("un ORF existe");
        assert_eq!(orf.start, 2);
        assert_eq!(orf.protein, "MA"); // le `*` est retiré
        assert_eq!(orf.end, 2 + 3 * 3); // AUG GCU UAA
    }

    #[test]
    fn find_orf_returns_none_without_stop() {
        assert!(find_orf("AUGGCUGCUGCU").is_none());
    }

    #[test]
    fn compare_simple_detects_the_sickle_cell_mutation() {
        // HBB : GAG (Glu) -> GTG (Val), position 5 du 2e codon.
        let normal = "ATGGTGCACCTGACTCCTGAGGAG";
        let sickle = "ATGGTGCACCTGACTCCTGTGGAG";
        let cmp = compare_simple(normal, sickle);
        assert_eq!(cmp.n_differences(), 1);
        assert_eq!(cmp.differences, vec![19]);
        assert!(cmp.identity_pct > 95.0);
    }

    #[test]
    fn compare_simple_handles_empty_input() {
        let cmp = compare_simple("", "");
        assert_eq!(cmp.length, 0);
        assert_eq!(cmp.identity_pct, 0.0);
    }
}
