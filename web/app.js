// DZ-Gene — منطق الواجهة (المرحلة 1)
// لا منطق علمي هنا: هذا الملف يعرض ويترجم ويربط النقرات فقط.
// كل علم بيولوجي يعيش في dzgene-core عبر أوامر Tauri.

const SUPPORTED = ["ar", "fr"];
const DEFAULT_LANG = "ar";

let translations = {};
let currentLang = DEFAULT_LANG;

// تحميل ملف ترجمة لغة ما (من مجلّد i18n).
async function loadLanguage(lang) {
  const res = await fetch(`i18n/${lang}.json`);
  if (!res.ok) throw new Error(`تعذّر تحميل لغة: ${lang}`);
  return res.json();
}

// تطبيق الترجمة على كل عنصر يحمل data-i18n.
function applyTranslations(dict) {
  document.querySelectorAll("[data-i18n]").forEach((el) => {
    const key = el.getAttribute("data-i18n");
    if (dict[key] !== undefined) el.textContent = dict[key];
  });

  // اتجاه الصفحة يتبع اللغة (rtl/ltr).
  document.documentElement.setAttribute("dir", dict.dir || "ltr");
  document.documentElement.setAttribute("lang", currentLang);

  // إبراز زرّ اللغة النشطة.
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.lang === currentLang);
  });
}

// تبديل اللغة فورًا.
async function switchLanguage(lang) {
  if (!SUPPORTED.includes(lang)) return;
  if (!translations[lang]) translations[lang] = await loadLanguage(lang);
  currentLang = lang;
  applyTranslations(translations[lang]);
}

// الاستجابة للنقر على عقدة في خريطة المفهوم.
function handleNodeClick(node) {
  if (node === "gene") {
    // مدخل المورثة: يعمل فعليًا (المرحلة 1). سنبنيه في الخطوة القادمة.
    console.log("فتح مدخل المورثة (التتابعات)");
    // مؤقتًا: رسالة حتى نبني الشاشة الفعلية.
    alert(translations[currentLang].node_gene + " — " + "قريبًا في الخطوة القادمة");
  } else {
    // البروتين والوظيفة: مُعلَّمان «قريبًا».
    const key = node === "protein" ? "node_protein" : "node_function";
    alert(translations[currentLang][key] + " — " + translations[currentLang].coming_soon);
  }
}

// ربط الأحداث بعد تحميل الصفحة.
function bindEvents() {
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.addEventListener("click", () => switchLanguage(btn.dataset.lang));
  });
  document.querySelectorAll(".node").forEach((node) => {
    node.addEventListener("click", () => handleNodeClick(node.dataset.node));
  });
}

// نقطة البداية.
async function main() {
  bindEvents();
  await switchLanguage(DEFAULT_LANG);
}

main();
