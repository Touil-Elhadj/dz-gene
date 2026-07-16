//! Distances évolutives — portage de `geniegen/bio/distances.py`.
//!
//! Trois modèles, du plus simple au plus corrigé. En 3AS on utilise surtout la
//! p-distance ; Jukes-Cantor et Kimura-2P servent à montrer aux élèves que
//! « compter les différences » sous-estime le nombre réel de substitutions
//! (mutations multiples au même site, ou retour à la base d'origine).

use serde::{Deserialize, Serialize};

/// Modèle de correction de distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Model {
    /// Proportion brute de sites différents.
    P,
    /// Jukes-Cantor : corrige les substitutions multiples (toutes équiprobables).
    JukesCantor,
    /// Kimura 2-paramètres : distingue transitions et transversions.
    Kimura2P,
}

/// Décompte des sites comparables entre deux séquences alignées.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SiteCounts {
    /// Sites comparés (gaps et `N` exclus).
    pub compared: usize,
    /// Transitions : A↔G (purines) ou C↔T (pyrimidines).
    pub transitions: usize,
    /// Transversions : purine ↔ pyrimidine.
    pub transversions: usize,
}

impl SiteCounts {
    /// Nombre total de sites qui diffèrent (transitions + transversions).
    pub fn differences(&self) -> usize {
        self.transitions + self.transversions
    }
}

#[inline]
fn is_purine(b: u8) -> bool {
    matches!(b, b'A' | b'G')
}

#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Compte transitions et transversions sur deux séquences **déjà alignées**.
///
/// Les positions contenant un gap ou une base ambiguë sont exclues du
/// dénominateur (« pairwise deletion »).
pub fn count_sites(a: &str, b: &str) -> SiteCounts {
    let mut c = SiteCounts::default();

    for (x, y) in a.bytes().zip(b.bytes()) {
        // `U` (ARN) est traité comme `T`.
        let x = if x.to_ascii_uppercase() == b'U' { b'T' } else { x.to_ascii_uppercase() };
        let y = if y.to_ascii_uppercase() == b'U' { b'T' } else { y.to_ascii_uppercase() };

        if !is_valid_base(x) || !is_valid_base(y) {
            continue;
        }
        c.compared += 1;

        if x == y {
            continue;
        }
        if is_purine(x) == is_purine(y) {
            c.transitions += 1;
        } else {
            c.transversions += 1;
        }
    }
    c
}

/// Distance évolutive entre deux séquences alignées, selon le modèle choisi.
///
/// Retourne `None` quand la distance n'est pas définie :
/// * aucun site comparable ;
/// * saturation — le modèle diverge (log d'un nombre ≤ 0). C'est un résultat
///   scientifique, pas une erreur : au-delà d'un certain niveau de divergence,
///   la correction n'a plus de sens. On préfère `None` à un `f64::INFINITY`
///   silencieux qui polluerait la matrice et l'arbre.
pub fn distance(a: &str, b: &str, model: Model) -> Option<f64> {
    let c = count_sites(a, b);
    if c.compared == 0 {
        return None;
    }
    let n = c.compared as f64;

    match model {
        Model::P => Some(c.differences() as f64 / n),

        Model::JukesCantor => {
            let p = c.differences() as f64 / n;
            let inner = 1.0 - (4.0 / 3.0) * p;
            (inner > 0.0).then(|| -0.75 * inner.ln())
        }

        Model::Kimura2P => {
            let s = c.transitions as f64 / n; // proportion de transitions
            let v = c.transversions as f64 / n; // proportion de transversions
            let t1 = 1.0 - 2.0 * s - v;
            let t2 = 1.0 - 2.0 * v;
            (t1 > 0.0 && t2 > 0.0).then(|| -0.5 * t1.ln() - 0.25 * t2.ln())
        }
    }
}

/// Matrice de distances symétrique, diagonale nulle.
///
/// Les distances indéfinies (saturation) sont représentées par `None`.
pub fn distance_matrix(seqs: &[String], model: Model) -> Vec<Vec<Option<f64>>> {
    let n = seqs.len();
    let mut m = vec![vec![Some(0.0); n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let d = distance(&seqs[i], &seqs[j], model);
            m[i][j] = d;
            m[j][i] = d;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }

    #[test]
    fn identical_sequences_have_zero_distance() {
        for m in [Model::P, Model::JukesCantor, Model::Kimura2P] {
            approx(distance("ATGC", "ATGC", m).unwrap(), 0.0);
        }
    }

    #[test]
    fn p_distance_is_the_raw_proportion() {
        // 1 différence sur 4 sites.
        approx(distance("ATGC", "ATGA", Model::P).unwrap(), 0.25);
    }

    #[test]
    fn transitions_and_transversions_are_distinguished() {
        // A->G : transition (purine->purine).
        let c = count_sites("AAAA", "GAAA");
        assert_eq!(c.transitions, 1);
        assert_eq!(c.transversions, 0);

        // A->C : transversion (purine->pyrimidine).
        let c = count_sites("AAAA", "CAAA");
        assert_eq!(c.transitions, 0);
        assert_eq!(c.transversions, 1);
    }

    #[test]
    fn gaps_are_excluded_from_the_denominator() {
        let c = count_sites("ATGC", "AT-C");
        assert_eq!(c.compared, 3); // le gap ne compte pas
        assert_eq!(c.differences(), 0);
    }

    #[test]
    fn jukes_cantor_exceeds_p_distance() {
        // La correction doit toujours *augmenter* la distance : c'est tout
        // l'intérêt pédagogique (les substitutions multiples sont invisibles).
        let a = "ATGCATGCATGCATGCATGC";
        let b = "AGGCTTGAATGCATCCTTGC";
        let p = distance(a, b, Model::P).unwrap();
        let jc = distance(a, b, Model::JukesCantor).unwrap();
        assert!(jc > p, "JC ({jc}) doit dépasser p ({p})");
    }

    #[test]
    fn saturation_yields_none_not_infinity() {
        // 100 % de différences : JC diverge (log d'un négatif).
        assert_eq!(distance("AAAA", "GGGG", Model::JukesCantor), None);
    }

    #[test]
    fn no_comparable_site_yields_none() {
        assert_eq!(distance("----", "----", Model::P), None);
        assert_eq!(distance("", "", Model::P), None);
    }

    #[test]
    fn rna_u_is_treated_as_t() {
        approx(distance("AUGC", "ATGC", Model::P).unwrap(), 0.0);
    }

    #[test]
    fn matrix_is_symmetric_with_null_diagonal() {
        let seqs: Vec<String> = ["ATGC", "ATGA", "TTGA"].iter().map(|s| s.to_string()).collect();
        let m = distance_matrix(&seqs, Model::P);

        for i in 0..3 {
            approx(m[i][i].unwrap(), 0.0);
            for j in 0..3 {
                assert_eq!(m[i][j], m[j][i]);
            }
        }
    }
}
