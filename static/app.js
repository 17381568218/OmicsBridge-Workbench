// Tab switching
function switchTab(id) {
  document.querySelectorAll('.tab').forEach(b => b.classList.remove('active'));
  const tabIds = ['tab1', 'tab6', 'tab2', 'tab3', 'tab4', 'tab5', 'tab7'];
  const idx = tabIds.indexOf(id);
  if (idx >= 0) document.querySelectorAll('.tab')[idx].classList.add('active');
  document.querySelectorAll('.tab-content').forEach(d => d.classList.remove('active'));
  document.getElementById(id).classList.add('active');
}

// Tab 1: File name display
['pos_com','pos_iden','neg_com','neg_iden'].forEach(id => {
  document.getElementById(id).addEventListener('change', e => {
    document.getElementById(id+'_name').textContent = e.target.files[0]?.name || '';
  });
});

document.getElementById('met_matrix')?.addEventListener('change', e => {
  document.getElementById('met_matrix_name').textContent = e.target.files[0]?.name || '';
});

// Tab 1: Load filter checkboxes
async function loadFilters() {
  try {
    const r = await fetch('/api/origins_tissues');
    const d = await r.json();
    buildCheckboxes('origin_chk', 'Keep Origins', d.origins, true);
    buildCheckboxes('tissue_chk', 'Priority Tissues', d.tissues, false);
  } catch(e) { console.error(e); }
}

function buildCheckboxes(id, label, values, defaultChecked) {
  const container = document.getElementById('filter_checkboxes');
  const div = document.createElement('div');
  div.className = 'chk-area';
  div.innerHTML = '<b style="font-size:0.85em">'+label+' (click to toggle):</b><br>'+
    '<div class="chk-actions"><button type="button" onclick="selAll(\''+id+'\')">All</button><button type="button" onclick="selNone(\''+id+'\')">None</button></div>'+
    '<div class="chk-group" id="'+id+'"></div>';
  container.appendChild(div);

  const group = document.getElementById(id);
  values.forEach(v => {
    const span = document.createElement('span');
    span.className = 'chk-tag' + (defaultChecked ? ' checked' : '');
    span.textContent = v;
    span.addEventListener('click', function() { this.classList.toggle('checked'); });
    group.appendChild(span);
  });
}

function selAll(id) {
  document.getElementById(id).querySelectorAll('.chk-tag').forEach(s => s.classList.add('checked'));
}
function selNone(id) {
  document.getElementById(id).querySelectorAll('.chk-tag').forEach(s => s.classList.remove('checked'));
}

function getChecked(id) {
  return Array.from(document.getElementById(id).querySelectorAll('.chk-tag.checked')).map(s => s.textContent.trim());
}

// Tab 1: Pipeline
let pollTimer = null;

async function runPipeline() {
  const posCom = document.getElementById('pos_com').files[0];
  const posIden = document.getElementById('pos_iden').files[0];
  if (!posCom || !posIden) { alert('POS COM and IDEN are required'); return; }

  const form = new FormData();
  form.append('pos_com', posCom);
  form.append('pos_iden', posIden);
  const negCom = document.getElementById('neg_com').files[0];
  const negIden = document.getElementById('neg_iden').files[0];
  if (negCom) form.append('neg_com', negCom);
  if (negIden) form.append('neg_iden', negIden);

  form.append('score_cutoff', document.getElementById('score_cut').value || '');
  form.append('frag_cutoff', document.getElementById('frag_cut').value || '');
  form.append('dccs_cutoff', document.getElementById('dccs_cut').value || '');
  form.append('allowed_origins', getChecked('origin_chk').join(','));
  form.append('priority_tissues', getChecked('tissue_chk').join(','));

  document.getElementById('btn_run').disabled = true;
  document.getElementById('btn_stop').disabled = false;
  document.getElementById('log_output').textContent = '';
  document.getElementById('downloads_section').style.display = 'none';

  try {
    const r = await fetch('/api/run', { method: 'POST', body: form });
    const d = await r.json();
    if (d.error) { alert(d.error); resetBtns(); return; }
    pollTask(d.task_id);
  } catch(e) { alert('Error: '+e); resetBtns(); }
}

function stopPipeline() {
  clearInterval(pollTimer);
  resetBtns();
}

function resetBtns() {
  document.getElementById('btn_run').disabled = false;
  document.getElementById('btn_stop').disabled = true;
}

function pollTask(taskId) {
  pollTimer = setInterval(async () => {
    try {
      const r = await fetch('/api/status/' + taskId);
      const d = await r.json();
      document.getElementById('log_output').textContent = d.log || '';
      document.getElementById('log_output').scrollTop = document.getElementById('log_output').scrollHeight;

      if (!d.running) {
        clearInterval(pollTimer);
        resetBtns();
        if (d.files && d.files.length > 0) {
          document.getElementById('downloads_section').style.display = 'block';
          const links = document.getElementById('download_links');
          links.innerHTML = d.files.map(f =>
            '<a class="dl-link" href="/api/download/'+taskId+'/'+f+'" download>📥 '+f+'</a>'
          ).join('');
        }
      }
    } catch(e) { clearInterval(pollTimer); resetBtns(); }
  }, 500);
}

// Tabs 2-5: Query
async function qRun(tabId, endpoint) {
  const ta = document.getElementById('ta_'+tabId);
  const fileInput = document.getElementById('file_'+tabId);
  const text = ta.value.trim();
  const file = fileInput.files[0];

  if (!text && !file) { alert('Paste text or select a file'); return; }

  const form = new FormData();
  if (file) {
    form.append('file', file);
  } else {
    form.append('text', text);
  }

  try {
    const r = await fetch(endpoint, { method: 'POST', body: form });
    const d = await r.json();
    if (d.error) { alert(d.error); return; }

    const tbl = document.getElementById('tbl_'+tabId);
    tbl.querySelector('thead').innerHTML = '<tr>'+d.columns.map(c => '<th>'+c+'</th>').join('')+'</tr>';
    tbl.querySelector('tbody').innerHTML = d.rows.map(row =>
      '<tr>'+row.map(v => '<td>'+(v||'')+'</td>').join('')+'</tr>'
    ).join('');
    document.getElementById('count_'+tabId).textContent = d.rows.length + ' results';

    // Store for export
    tbl._data = d;
  } catch(e) { alert('Error: '+e); }
}

function qImport(tabId) {
  const file = document.getElementById('file_'+tabId).files[0];
  if (!file) return;
  // Auto-trigger query on file select
  qRun(tabId, QUERY_TABS.find(t => t.id === tabId).endpoint);
}

function qExport(tabId, colsStr) {
  const tbl = document.getElementById('tbl_'+tabId);
  const data = tbl._data;
  if (!data || !data.rows.length) { alert('No results to export'); return; }

  let csv = data.columns.join(',') + '\n';
  csv += data.rows.map(row => row.map(v => '"'+(v||'').replace(/"/g,'""')+'"').join(',')).join('\n');

  const blob = new Blob(['﻿'+csv], { type: 'text/csv;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'export.csv'; a.click();
  URL.revokeObjectURL(url);
}

// Tab 6: Metabolomics analysis
let metPollTimer = null;
const MET_COLORS = ['#2563eb', '#dc2626', '#059669', '#d97706', '#7c3aed', '#0891b2', '#be123c', '#4d7c0f'];

async function runMetabolomics() {
  const matrix = document.getElementById('met_matrix').files[0];
  if (!matrix) { alert('Please select a metabolomics matrix CSV'); return; }

  const form = new FormData();
  form.append('matrix', matrix);
  form.append('group_a', document.getElementById('met_group_a').value.trim());
  form.append('group_b', document.getElementById('met_group_b').value.trim());

  document.getElementById('btn_met_run').disabled = true;
  document.getElementById('met_log_output').textContent = '';
  document.getElementById('met_downloads_section').style.display = 'none';
  document.getElementById('met_metrics').textContent = 'Running analysis...';

  try {
    const r = await fetch('/api/metabolomics/run', { method: 'POST', body: form });
    const d = await r.json();
    if (d.error) { alert(d.error); document.getElementById('btn_met_run').disabled = false; return; }
    pollMetabolomics(d.task_id);
  } catch (e) {
    alert('Error: ' + e);
    document.getElementById('btn_met_run').disabled = false;
  }
}

function pollMetabolomics(taskId) {
  clearInterval(metPollTimer);
  metPollTimer = setInterval(async () => {
    try {
      const r = await fetch('/api/metabolomics/status/' + taskId);
      const d = await r.json();
      document.getElementById('met_log_output').textContent = d.log || '';
      document.getElementById('met_log_output').scrollTop = document.getElementById('met_log_output').scrollHeight;

      if (!d.running) {
        clearInterval(metPollTimer);
        document.getElementById('btn_met_run').disabled = false;
        if (d.error) { document.getElementById('met_metrics').textContent = d.error; return; }
        if (d.result) renderMetabolomics(d.result, taskId, d.files || []);
      }
    } catch (e) {
      clearInterval(metPollTimer);
      document.getElementById('btn_met_run').disabled = false;
    }
  }, 700);
}

function renderMetabolomics(result, taskId, files) {
  renderMetMetrics(result);
  renderSampleTable(result.summary.samples || []);
  renderPca(result.pca || { points: [] });
  renderVolcano(result.volcano || []);
  renderHeatmap(result.heatmap || {});
  renderTopTable(result.top_table || []);

  const section = document.getElementById('met_downloads_section');
  const links = document.getElementById('met_download_links');
  section.style.display = files.length ? 'block' : 'none';
  links.innerHTML = files.map(f =>
    '<a class="dl-link" href="/api/download/' + taskId + '/' + f + '" download>' + f + '</a>'
  ).join('');
}

function renderMetMetrics(result) {
  const s = result.summary;
  const qc30 = result.qc.features_with_qc_rsd_le_30;
  const html = [
    ['Features', s.feature_count],
    ['Samples', s.sample_count],
    ['Biological', s.biological_sample_count],
    ['QC', s.qc_sample_count],
    ['Zero %', s.zero_percent],
    ['Missing %', s.missing_percent],
    ['Significant', result.significant_count],
    ['QC RSD <= 30%', qc30 === null ? 'NA' : qc30],
  ].map(([k, v]) => '<div class="metric"><span>' + k + '</span><strong>' + v + '</strong></div>').join('');
  document.getElementById('met_metrics').classList.remove('empty-state');
  document.getElementById('met_metrics').innerHTML = html;

  document.getElementById('met_warnings').innerHTML = (s.warnings || []).map(w =>
    '<div class="warning-item">' + escapeHtml(w) + '</div>'
  ).join('');
}

function renderSampleTable(samples) {
  const cols = ['sample_id', 'mode', 'group', 'replicate', 'is_qc'];
  const tbl = document.getElementById('met_sample_table');
  tbl.querySelector('thead').innerHTML = '<tr>' + cols.map(c => '<th>' + c + '</th>').join('') + '</tr>';
  tbl.querySelector('tbody').innerHTML = samples.map(row =>
    '<tr>' + cols.map(c => '<td>' + escapeHtml(String(row[c] ?? '')) + '</td>').join('') + '</tr>'
  ).join('');
}

function renderTopTable(rows) {
  const cols = ['Compound_ID', 'QC_RSD_percent', 'log2FC', 'p_value', 'FDR', 'direction'];
  const tbl = document.getElementById('met_top_table');
  tbl.querySelector('thead').innerHTML = '<tr>' + cols.map(c => '<th>' + c + '</th>').join('') + '</tr>';
  tbl.querySelector('tbody').innerHTML = rows.map(row =>
    '<tr>' + cols.map(c => '<td>' + formatCell(row[c]) + '</td>').join('') + '</tr>'
  ).join('');
}

function renderPca(pca) {
  const points = pca.points || [];
  if (!points.length) { document.getElementById('met_pca').innerHTML = '<div class="empty-state">No PCA data</div>'; return; }
  const xs = points.map(p => p.pc1), ys = points.map(p => p.pc2);
  const scale = makeScale(xs, ys, 42, 318, 222, 32);
  const groups = Array.from(new Set(points.map(p => p.group)));
  const color = g => MET_COLORS[groups.indexOf(g) % MET_COLORS.length];
  let svg = chartSvgStart();
  svg += axisSvg('PC1 ' + (pca.explained?.[0] || 0) + '%', 'PC2 ' + (pca.explained?.[1] || 0) + '%');
  points.forEach(p => {
    svg += '<circle cx="' + scale.x(p.pc1) + '" cy="' + scale.y(p.pc2) + '" r="6" fill="' + color(p.group) + '"><title>' + escapeHtml(p.sample) + '</title></circle>';
  });
  svg += legendSvg(groups, color);
  svg += '</svg>';
  document.getElementById('met_pca').innerHTML = svg;
}

function renderVolcano(rows) {
  if (!rows.length) { document.getElementById('met_volcano').innerHTML = '<div class="empty-state">No volcano data</div>'; return; }
  const xs = rows.map(r => r.log2FC), ys = rows.map(r => r.negLog10P);
  const scale = makeScale(xs, ys, 42, 318, 222, 32);
  let svg = chartSvgStart();
  svg += axisSvg('log2FC', '-log10(P)');
  rows.forEach(r => {
    const fill = r.significant ? (r.log2FC > 0 ? '#dc2626' : '#2563eb') : '#94a3b8';
    const radius = r.significant ? 4 : 2.5;
    svg += '<circle cx="' + scale.x(r.log2FC) + '" cy="' + scale.y(r.negLog10P) + '" r="' + radius + '" fill="' + fill + '" opacity="0.78"><title>' + escapeHtml(r.id) + '</title></circle>';
  });
  svg += '</svg>';
  document.getElementById('met_volcano').innerHTML = svg;
}

function renderHeatmap(h) {
  if (!h.features?.length || !h.samples?.length) {
    document.getElementById('met_heatmap').innerHTML = '<div class="empty-state">No heatmap data</div>';
    return;
  }
  const cellW = Math.max(18, Math.min(34, Math.floor(620 / h.samples.length)));
  const cellH = 17;
  const html = h.features.map((feature, i) => {
    const cells = h.matrix[i].map(v => '<span class="hm-cell" style="width:' + cellW + 'px;height:' + cellH + 'px;background:' + heatColor(v) + '" title="' + feature + ': ' + v + '"></span>').join('');
    return '<div class="hm-row"><span class="hm-label">' + escapeHtml(feature) + '</span><span class="hm-cells">' + cells + '</span></div>';
  }).join('');
  const sampleLabels = '<div class="hm-samples"><span class="hm-label"></span>' + h.samples.map(s => '<span style="width:' + cellW + 'px">' + escapeHtml(shortSample(s)) + '</span>').join('') + '</div>';
  document.getElementById('met_heatmap').innerHTML = '<div class="hm-scroll">' + sampleLabels + html + '</div>';
}

function chartSvgStart() {
  return '<svg viewBox="0 0 360 260" class="svg-chart" role="img">';
}

function axisSvg(xLabel, yLabel) {
  return '<line x1="42" y1="222" x2="326" y2="222" stroke="#cbd5e1"/>' +
    '<line x1="42" y1="222" x2="42" y2="26" stroke="#cbd5e1"/>' +
    '<text x="184" y="250" text-anchor="middle" class="axis-label">' + escapeHtml(xLabel) + '</text>' +
    '<text x="14" y="124" text-anchor="middle" transform="rotate(-90 14 124)" class="axis-label">' + escapeHtml(yLabel) + '</text>';
}

function legendSvg(groups, color) {
  return groups.map((g, i) => {
    const y = 20 + i * 16;
    return '<circle cx="282" cy="' + y + '" r="4" fill="' + color(g) + '"/><text x="291" y="' + (y + 4) + '" class="legend-label">' + escapeHtml(g) + '</text>';
  }).join('');
}

function makeScale(xs, ys, x0, x1, y0, y1) {
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const padX = (maxX - minX || 1) * 0.1;
  const padY = (maxY - minY || 1) * 0.1;
  return {
    x: v => x0 + ((v - minX + padX) / (maxX - minX + padX * 2)) * (x1 - x0),
    y: v => y0 - ((v - minY + padY) / (maxY - minY + padY * 2)) * (y0 - y1),
  };
}

function heatColor(v) {
  const x = Math.max(-3, Math.min(3, Number(v) || 0)) / 3;
  if (x >= 0) {
    const b = Math.round(245 - x * 150), g = Math.round(245 - x * 95), r = 220;
    return 'rgb(' + r + ',' + g + ',' + b + ')';
  }
  const a = Math.abs(x), r = Math.round(245 - a * 150), g = Math.round(245 - a * 75), b = 210;
  return 'rgb(' + r + ',' + g + ',' + b + ')';
}

function formatCell(v) {
  if (v === null || v === undefined) return '';
  if (typeof v === 'number') {
    if (Math.abs(v) < 0.001 && v !== 0) return v.toExponential(2);
    return String(Math.round(v * 10000) / 10000);
  }
  return escapeHtml(String(v));
}

function shortSample(s) {
  const parts = String(s).split('-');
  return parts.length >= 2 ? parts.slice(-2).join('-') : String(s);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

// Init
loadFilters();
