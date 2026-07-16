//! Couche de commandes Tauri — **volontairement rachitique**.
//!
//! Aucune logique scientifique ici. Chaque commande se contente de :
//! 1. recevoir des données du front (chaînes, entiers) ;
//! 2. appeler `dzgene_core` ;
//! 3. renvoyer une structure sérialisable.
//!
//! Si vous vous surprenez à écrire un `if` biologique dans ce fichier, la
//! fonction est au mauvais endroit : elle appartient à `dzgene-core`, où elle
//! sera testable sans lancer d'application.

use dzgene_core::{distances, edi, genetics};
use serde::Serialize;

/// Erreur remontée au front-end.
#[derive(Debug, Serialize)]
pub struct CmdError {
    /// Clé i18n (le front choisit la langue : FR / AR / EN).
    pub key: String,
    /// Détail technique, pour le journal — pas pour l'élève.
    pub detail: String,
}

impl CmdError {
    fn new(key: &str, detail: impl std::fmt::Display) -> Self {
        Self { key: key.into(), detail: detail.to_string() }
    }
}

type CmdResult<T> = Result<T, CmdError>;

// --------------------------------------------------------------------------
// Génétique
// --------------------------------------------------------------------------

#[tauri::command]
pub fn transcribe(template_dna: &str) -> String {
    genetics::transcribe(template_dna)
}

#[tauri::command]
pub fn reverse_complement(dna: &str) -> String {
    genetics::reverse_complement(dna)
}

#[tauri::command]
pub fn translate(seq: &str, frame: usize, to_stop: bool) -> CmdResult<genetics::Translation> {
    if frame > 2 {
        return Err(CmdError::new("error.invalid_frame", frame));
    }
    Ok(genetics::translate(seq, frame, to_stop))
}

#[tauri::command]
pub fn find_orf(seq: &str) -> Option<genetics::Orf> {
    genetics::find_orf(seq)
}

#[tauri::command]
pub fn compare_simple(a: &str, b: &str) -> genetics::SimpleComparison {
    genetics::compare_simple(a, b)
}

// --------------------------------------------------------------------------
// Banque .edi
// --------------------------------------------------------------------------

/// Ouvre un fichier `.edi` depuis le disque.
///
/// C'est ici — et **pas** dans `dzgene-core` — que vit l'accès au système de
/// fichiers. La crate scientifique reste pure et testable hors-ligne.
#[tauri::command]
pub fn open_edi(path: &str) -> CmdResult<edi::EdiFile> {
    let bytes = std::fs::read(path).map_err(|e| CmdError::new("error.read_failed", e))?;
    let file = edi::parse_bytes(&bytes);

    if file.is_empty() {
        return Err(CmdError::new("error.empty_bank", path));
    }
    Ok(file)
}

#[tauri::command]
pub fn save_edi(path: &str, file: edi::EdiFile) -> CmdResult<()> {
    std::fs::write(path, edi::to_bytes(&file))
        .map_err(|e| CmdError::new("error.write_failed", e))
}

// --------------------------------------------------------------------------
// Phylogénie
// --------------------------------------------------------------------------

/// Matrice de distances.
///
/// `None` (→ `null` en JSON) signale une distance **indéfinie par saturation** :
/// le front doit l'afficher comme telle (« ∞ » / case grisée) et non comme un
/// zéro. C'est un point pédagogique, pas un détail d'implémentation.
#[tauri::command]
pub fn distance_matrix(seqs: Vec<String>, model: distances::Model) -> Vec<Vec<Option<f64>>> {
    distances::distance_matrix(&seqs, model)
}

#[tauri::command]
pub fn version() -> &'static str {
    dzgene_core::VERSION
}
