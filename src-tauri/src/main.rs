//! DZ-Gene — point d'entrée de l'application de bureau (Tauri).
//!
//! Ce binaire ne fait qu'assembler : il enregistre les commandes et lance la
//! webview. Toute la science est dans `dzgene-core` ; toute la présentation est
//! dans `web/`. Ce fichier ne doit jamais grossir.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod commands;

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            commands::transcribe,
            commands::reverse_complement,
            commands::translate,
            commands::find_orf,
            commands::compare_simple,
            commands::open_edi,
            commands::save_edi,
            commands::distance_matrix,
            commands::version,
        ])
        .run(tauri::generate_context!())
        .expect("échec du démarrage de DZ-Gene");
}
