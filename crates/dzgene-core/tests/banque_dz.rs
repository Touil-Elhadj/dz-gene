//! Test de **non-régression** contre la banque DZ réelle.
//!
//! C'est le filet de sécurité de toute la migration : si ces tests passent, un
//! enseignant peut ouvrir dans la version Rust exactement les mêmes fichiers
//! qu'il ouvrait dans la version PySide6. Aucune fonctionnalité brillante ne
//! rachète une banque qui ne s'ouvre plus.
//!
//! Les fixtures sont des copies **non modifiées** de `banque-dz/`, couvrant les
//! programmes 1AS, 2AS et 3AS.

use dzgene_core::edi::{self, SeqType};
use dzgene_core::genetics;
use std::fs;
use std::path::PathBuf;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn all_fixtures() -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(fixtures_dir())
        .expect("le dossier de fixtures doit exister")
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "edi"))
        .collect();
    files.sort();
    assert!(!files.is_empty(), "aucune fixture .edi trouvée");
    files
}

#[test]
fn every_bank_file_parses_without_panicking() {
    for path in all_fixtures() {
        let bytes = fs::read(&path).unwrap();
        let file = edi::parse_bytes(&bytes);
        assert!(
            !file.is_empty(),
            "{} : aucune séquence lue",
            path.display()
        );
    }
}

#[test]
fn every_sequence_has_a_name_and_residues() {
    for path in all_fixtures() {
        let file = edi::parse_bytes(&fs::read(&path).unwrap());
        for seq in file.iter() {
            assert!(
                !seq.name.trim().is_empty(),
                "{} : séquence sans nom",
                path.display()
            );
            assert!(
                seq.length() > 0,
                "{} : séquence « {} » vide",
                path.display(),
                seq.name
            );
        }
    }
}

#[test]
fn proteins_are_normalised_to_one_letter_code() {
    // Toute protéine, quel que soit son `Type` d'origine (2 ou 3), doit sortir
    // du parseur en code 1 lettre. C'est l'invariant sur lequel s'appuient
    // l'alignement et l'affichage.
    for path in all_fixtures() {
        let file = edi::parse_bytes(&fs::read(&path).unwrap());
        for seq in file.iter().filter(|s| s.seqtype.is_protein()) {
            for c in seq.residues.chars() {
                assert!(
                    c.is_ascii_uppercase() || c == '*' || c == '-',
                    "{} / {} : résidu inattendu « {c} » (protéine non normalisée ?)",
                    path.display(),
                    seq.name
                );
            }
        }
    }
}

#[test]
fn nucleic_sequences_contain_only_valid_symbols() {
    for path in all_fixtures() {
        let file = edi::parse_bytes(&fs::read(&path).unwrap());
        for seq in file.iter().filter(|s| s.seqtype == SeqType::Adn) {
            for c in seq.residues.chars() {
                assert!(
                    c.is_ascii_uppercase() || c == '-' || c == '*',
                    "{} / {} : symbole nucléique invalide « {c} »",
                    path.display(),
                    seq.name
                );
            }
        }
    }
}

#[test]
fn round_trip_is_stable_on_the_whole_bank() {
    // parse -> write -> parse doit être un point fixe. Sans cette garantie,
    // un enseignant qui ré-enregistre une séquence la corrompt silencieusement.
    for path in all_fixtures() {
        let original = edi::parse_bytes(&fs::read(&path).unwrap());
        let rewritten = edi::parse_bytes(&edi::to_bytes(&original));

        assert_eq!(
            original.len(),
            rewritten.len(),
            "{} : nombre de séquences modifié par l'aller-retour",
            path.display()
        );

        for (a, b) in original.iter().zip(rewritten.iter()) {
            assert_eq!(
                a.residues,
                b.residues,
                "{} / {} : résidus altérés par l'aller-retour",
                path.display(),
                a.name
            );
            assert_eq!(a.name, b.name, "{} : nom altéré", path.display());
            assert_eq!(a.seqtype, b.seqtype, "{} : type altéré", path.display());
        }
    }
}

#[test]
fn beta_globin_dna_translates_to_a_plausible_protein() {
    // Test scientifique de bout en bout : on lit un vrai fichier de la banque
    // 2AS, on traduit, et on vérifie que la protéine commence par une
    // méthionine — ce que tout élève doit retrouver à la main.
    let path = fixtures_dir().join("2AS - Expression genetique (ADN-ARN-Proteine)__GlobineBetaADN.edi");
    let file = edi::parse_bytes(&fs::read(&path).expect("fixture GlobineBetaADN présente"));

    let seq = file.iter().next().expect("au moins une séquence");
    let tr = genetics::translate(&seq.residues, 0, true);

    assert!(
        tr.protein.starts_with('M'),
        "la bêta-globine doit démarrer par une méthionine, obtenu : {}",
        &tr.protein.chars().take(10).collect::<String>()
    );
}
