// DZ-Gene — منطق الواجهة (المرحلة 1)
// لا منطق علمي هنا: العرض والترجمة والربط فقط.
// كل علم بيولوجي يعيش في dzgene-core عبر أوامر Tauri.

const SUPPORTED = ["ar", "fr"];
const DEFAULT_LANG = "ar";
const RESIDUES_PER_LINE = 60; // مثل Anagène تقريبًا

let translations = {};
let currentLang = DEFAULT_LANG;
let currentFile = null;   // EdiFile المفتوح
let currentSeqIndex = -1; // التتابع المختار

// ---- جسر الاستدعاء إلى Rust ----
// نستعمل window.__TAURI__ (بلا حزم npm) ليبقى الكود هجينًا.
async function callRust(cmd, args) {
  if (window.__TAURI__ && window.__TAURI__.core) {
    return window.__TAURI__.core.invoke(cmd, args);
  }
  // بديل المتصفح (لاحقًا): حاليًّا ننبّه أننا خارج Tauri.
  throw new Error("خارج بيئة Tauri: الاستدعاء غير متاح بعد");
}

// ---- الترجمة ----
async function loadLanguage(lang) {
  const res = await fetch(`i18n/${lang}.json`);
  if (!res.ok) throw new Error(`تعذّر تحميل لغة: ${lang}`);
  return res.json();
}

function t(key) {
  return (translations[currentLang] && translations[currentLang][key]) || key;
}

function applyTranslations() {
  const dict = translations[currentLang];
  document.querySelectorAll("[data-i18n]").forEach((el) => {
    const key = el.getAttribute("data-i18n");
    if (dict[key] !== undefined) el.textContent = dict[key];
  });
  document.documentElement.setAttribute("dir", dict.dir || "ltr");
  document.documentElement.setAttribute("lang", currentLang);
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.lang === currentLang);
  });
}

async function switchLanguage(lang) {
  if (!SUPPORTED.includes(lang)) return;
  if (!translations[lang]) translations[lang] = await loadLanguage(lang);
  currentLang = lang;
  applyTranslations();
  // إعادة رسم التتابع المعروض (المسطرة ثابتة، لكن النصوص المحيطة تتغيّر).
  if (currentSeqIndex >= 0) renderSequence(currentSeqIndex);
}

// ---- التنقّل بين الشاشات ----
function showScreen(name) {
  const map = document.querySelector(".concept-map");
  const geneScreen = document.querySelector(".screen-gene");
  if (name === "gene") {
    map.classList.add("is-hidden");
    geneScreen.hidden = false;
  } else {
    map.classList.remove("is-hidden");
    geneScreen.hidden = true;
  }
}

function handleNodeClick(node) {
  if (node === "gene") {
    showScreen("gene");
  } else {
    const key = node === "protein" ? "node_protein" : "node_function";
    alert(t(key) + " — " + t("coming_soon"));
  }
}

// ---- فتح ملف .edi عبر نواة Rust ----
async function openEdiFile() {
  try {
    // نافذة اختيار ملف (إضافة dialog).
    const path = await window.__TAURI__.dialog.open({
      multiple: false,
      filters: [{ name: "Anagène EDI", extensions: ["edi"] }],
    });
    if (!path) return; // ألغى المستخدم

    // نستدعي أمر Rust: open_edi يقرأ ويحلّل الملف.
    currentFile = await callRust("open_edi", { path });
    currentSeqIndex = -1;

    // اسم الملف في الشريط العلوي.
    const name = String(path).split(/[\\/]/).pop();
    document.getElementById("fileName").textContent = name;

    renderSeqList();
  } catch (err) {
    alert("خطأ في فتح الملف:\n" + (err.message || err));
  }
}

// ---- قائمة التتابعات على الجانب ----
function renderSeqList() {
  const list = document.getElementById("seqList");
  list.innerHTML = "";

  currentFile.sequences.forEach((seq, i) => {
    const item = document.createElement("div");
    item.className = "seq-item";
    item.innerHTML =
      `<div>${escapeHtml(seq.name)}</div>` +
      `<div class="seq-item-meta">${seq.residues.length} — ${seq.seqtype}</div>`;
    item.addEventListener("click", () => renderSequence(i));
    list.appendChild(item);
  });
}

// ---- عرض تتابع مع مسطرة الترقيم ----
function renderSequence(index) {
  currentSeqIndex = index;
  const seq = currentFile.sequences[index];

  // إبراز المختار في القائمة.
  document.querySelectorAll(".seq-item").forEach((el, i) => {
    el.classList.toggle("active", i === index);
  });

  const detail = document.getElementById("seqDetail");
  const isNucleic = seq.seqtype === "adn";

  detail.innerHTML =
    `<div class="seq-header">` +
      `<span class="seq-title">${escapeHtml(seq.name)}</span>` +
      `<span class="seq-meta">${t("length")}: ${seq.residues.length}</span>` +
      `<span class="seq-meta">${t("type_label")}: ${seq.seqtype}</span>` +
    `</div>` +
    `<div class="seq-block">${buildRuledSequence(seq.residues)}</div>` +
    (isNucleic
      ? `<div class="translate-row">` +
          `<button class="translate-btn" data-action="translate">${t("translate_btn")}</button>` +
          `<div class="protein-out" id="proteinOut" hidden></div>` +
        `</div>`
      : "");
}

// يبني المسطرة + التتابع بأسطر، بمحاذاة أحادية المسافة.
function buildRuledSequence(residues) {
  let html = "";
  for (let start = 0; start < residues.length; start += RESIDUES_PER_LINE) {
    const chunk = residues.slice(start, start + RESIDUES_PER_LINE);
    html += `<div class="ruler-line">${buildRuler(start, chunk.length)}</div>`;
    html += `<div class="residues-line">${chunk}</div>`;
  }
  return html;
}

// مسطرة: رقم كل 10 (يبدأ العدّ من 1)، وعلامة '.' عند كل عاشرة.
function buildRuler(offset, len) {
  let ruler = "";
  for (let i = 0; i < len; i++) {
    const pos = offset + i + 1; // 1-based (اصطلاح البيولوجيا)
    if (pos % 10 === 0) {
      const label = String(pos);
      // نُرجع المؤشّر إلى الوراء لنُحاذي الرقم بحيث ينتهي عند العاشرة.
      ruler = ruler.slice(0, ruler.length - (label.length - 1)) + label;
    } else {
      ruler += " ";
    }
  }
  return ruler;
}

function escapeHtml(s) {
  const div = document.createElement("div");
  div.textContent = s;
  return div.innerHTML;
}

// ---- الترجمة عبر نواة Rust ----
async function translateCurrent() {
  const seq = currentFile.sequences[currentSeqIndex];
  try {
    const result = await callRust("translate", {
      seq: seq.residues,
      frame: 0,
      toStop: false,
    });
    const out = document.getElementById("proteinOut");
    out.hidden = false;
    out.textContent = result.protein;
  } catch (err) {
    alert("خطأ في الترجمة:\n" + (err.message || err));
  }
}

// ---- ربط الأحداث ----
function bindEvents() {
  document.querySelectorAll(".lang-btn").forEach((btn) => {
    btn.addEventListener("click", () => switchLanguage(btn.dataset.lang));
  });
  document.querySelectorAll(".node").forEach((node) => {
    node.addEventListener("click", () => handleNodeClick(node.dataset.node));
  });

  // أحداث مفوّضة (delegation) لأن بعض الأزرار تُنشأ ديناميكيًّا.
  document.addEventListener("click", (e) => {
    const action = e.target.closest("[data-action]")?.dataset.action;
    if (action === "back") showScreen("map");
    else if (action === "open-edi") openEdiFile();
    else if (action === "import-global") {
      // ثابت offline-first: لا اتصال هنا. الميزة تأتي في مرحلتها.
      alert(t("import_global") + " — " + t("coming_soon"));
    }
    else if (action === "translate") translateCurrent();
  });
}

async function main() {
  bindEvents();
  await switchLanguage(DEFAULT_LANG);
}

main();
