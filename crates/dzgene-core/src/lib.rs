//! # dzgene-core
//!
//! Moteur scientifique **pur** de DZ-Gene — outil de génétique moléculaire
//! destiné aux élèves de 3AS (Sciences expérimentales) en Algérie.
//!
//! ## Contrat architectural
//!
//! Cette crate ne connaît **ni Tauri, ni le web, ni aucune interface**. Elle
//! n'ouvre aucun fichier et n'accède pas au réseau : elle reçoit des chaînes,
//! rend des structures. C'est ce qui la rend :
//!
//! * **testable** sans écran ni environnement graphique ;
//! * **réutilisable** — la même crate alimente l'application de bureau (Tauri)
//!   et la version navigateur (WebAssembly) ;
//! * **citable** — publiable seule sur crates.io et archivable sur Zenodo.
//!
//! La règle de dépendance est stricte et à sens unique :
//!
//! ```text
//! interface (Tauri / web)  →  dzgene-core
//! ```
//!
//! Jamais l'inverse. Si une fonction d'ici a besoin de connaître l'interface,
//! c'est qu'elle n'est pas à sa place.
//!
//! ## Principe directeur
//!
//! DZ-Gene doit rester **utilisable hors-ligne, en salle de classe**, sans
//! binaire tiers ni connexion. Chaque opération pédagogique est donc
//! implémentée ici en Rust pur.
//!
//! ## Exemple
//!
//! ```
//! use dzgene_core::genetics;
//!
//! // La mutation de la drépanocytose : GAG (Glu) -> GTG (Val).
//! let normal = "ATGGTGCACCTGACTCCTGAGGAG";
//! let sickle = "ATGGTGCACCTGACTCCTGTGGAG";
//!
//! let cmp = genetics::compare_simple(normal, sickle);
//! assert_eq!(cmp.n_differences(), 1);
//!
//! let p1 = genetics::translate(normal, 0, false).protein;
//! let p2 = genetics::translate(sickle, 0, false).protein;
//! assert_eq!(p1, "MVHLTPEE");
//! assert_eq!(p2, "MVHLTPVE"); // un seul acide aminé change
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod distances;
pub mod edi;
pub mod genetics;

/// Version de la crate, exposée à l'interface (barre d'état, rapports).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
