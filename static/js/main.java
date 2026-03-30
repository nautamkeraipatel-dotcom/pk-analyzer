/**
 * main.js — PK Bioavailability Analyzer Frontend
 * Handles: second-route toggle, form validation, API call, result rendering
 */
(function () {
  'use strict';
  // ── DOM refs ──────────────────────────────────────────────────────────────                                   const togBtn    = document.getElementById('second_toggle');
  const cardR2    = document.getElementById('card-r2');
  const btnGen    = document.getElementById('btn_generate');
  const errBox    = document.getElementById('error_box');
  const resIdle   = document.getElementById('result_idle');
  const resLoad   = document.getElementById('result_loading');
  const resOut    = document.getElementById('result_output');
  const resImg    = document.getElementById('result_img');
  const resLabel  = document.getElementById('result_label');
  const btnDl     = document.getElementById('btn_download');

  let secondRouteOn = false;
  let lastFilename  = 'pk_analysis.jpg';
  let lastImgB64    = null;

  // ── Toggle second route ───────────────────────────────────────────────────
  togBtn.addEventListener('click', () => {
    secondRouteOn = !secondRouteOn;
    togBtn.setAttribute('aria-pressed', secondRouteOn.toString());
    togBtn.querySelector('.toggle-text').textContent = secondRouteOn ? 'ON' : 'OFF';
    cardR2.classList.toggle('visible', secondRouteOn);
    cardR2.setAttribute('aria-hidden', (!secondRouteOn).toString());
  });

  // ── Helpers ───────────────────────────────────────────────────────────────
  function showError(msg) {
    errBox.textContent = '⚠  ' + msg;
    errBox.classList.add('visible');
    errBox.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }

  function clearError() {
    errBox.textContent = '';
    errBox.classList.remove('visible');
  }

  function setState(state) {
    // state: 'idle' | 'loading' | 'result'
    resIdle.style.display  = state === 'idle'    ? '' : 'none';
    resLoad.style.display  = state === 'loading' ? '' : 'none';
    resOut.style.display   = state === 'result'  ? '' : 'none';
    resLoad.classList.toggle('active', state === 'loading');
    resOut.classList.toggle('active',  state === 'result');
  }

  function getVal(id) {
    return document.getElementById(id).value.trim();
  }

  function countNums(str) {
    return str.replace(/,/g, ' ').split(/\s+/).filter(s => s !== '' && !isNaN(Number(s))).length;
  }

  // ── Validation ────────────────────────────────────────────────────────────
  function validate() {
    const drugName   = getVal('drug_name');
    const mec        = getVal('mec');
    const msc        = getVal('msc');
    const route1     = getVal('route1_name');
    const dose1      = getVal('dose1');
    const time1      = getVal('time1_str');
    const conc1      = getVal('conc1_str');

    if (!drugName)              return 'Drug name is required.';
    if (!route1)                return 'Primary route name is required.';
    if (!dose1 || isNaN(Number(dose1)) || Number(dose1) <= 0)
                                return 'Primary route: enter a valid dose (mg).';
    if (!time1)                 return 'Primary route: time points are required.';
    if (!conc1)                 return 'Primary route: plasma concentrations are required.';

    const t1n = countNums(time1);
    const c1n = countNums(conc1);
    if (t1n < 2)                return 'Primary route: need at least 2 time points.';
    if (t1n !== c1n)            return `Primary route: ${t1n} time point(s) but ${c1n} concentration value(s). Counts must match.`;

    if (parseFloat(mec) > 0 && parseFloat(msc) > 0) {
      if (parseFloat(mec) >= parseFloat(msc))
        return 'MEC must be less than MSC.';
    }

    if (secondRouteOn) {
      const route2 = getVal('route2_name');
      const dose2  = getVal('dose2');
      const time2  = getVal('time2_str');
      const conc2  = getVal('conc2_str');

      if (!route2)              return 'Secondary route name is required.';
      if (!dose2 || isNaN(Number(dose2)) || Number(dose2) <= 0)
                                return 'Secondary route: enter a valid dose (mg).';
      if (!time2)               return 'Secondary route: time points are required.';
      if (!conc2)               return 'Secondary route: plasma concentrations are required.';

      const t2n = countNums(time2);
      const c2n = countNums(conc2);
      if (t2n < 2)              return 'Secondary route: need at least 2 time points.';
      if (t2n !== c2n)          return `Secondary route: ${t2n} time point(s) but ${c2n} concentration value(s). Counts must match.`;
    }

    return null; // valid
  }

  // ── Generate ──────────────────────────────────────────────────────────────
  btnGen.addEventListener('click', async () => {
    clearError();

    const validationError = validate();
    if (validationError) {
      showError(validationError);
      return;
    }

    const payload = {
      drug_name:        getVal('drug_name'),
      mec:              parseFloat(getVal('mec'))  || 0,
      msc:              parseFloat(getVal('msc'))  || 0,
      route1_name:      getVal('route1_name'),
      dose1:            parseFloat(getVal('dose1')),
      time1_str:        getVal('time1_str'),
      conc1_str:        getVal('conc1_str'),
      has_second_route: secondRouteOn,
    };

    if (secondRouteOn) {
      payload.route2_name = getVal('route2_name');
      payload.dose2       = parseFloat(getVal('dose2'));
      payload.time2_str   = getVal('time2_str');
      payload.conc2_str   = getVal('conc2_str');
    }

    btnGen.disabled = true;
    btnGen.querySelector('.btn-label').textContent = 'Analyzing…';
    setState('loading');

    try {
      const resp = await fetch('/analyze', {
        method:  'POST',
        headers: { 'Content-Type': 'application/json' },
        body:    JSON.stringify(payload),
      });

      const data = await resp.json();

      if (!resp.ok || data.error) {
        showError(data.error || 'Server returned an error. Check your inputs.');
        setState('idle');
        return;
      }

      // Render result
      lastImgB64   = data.image;
      lastFilename = data.filename;

      resImg.src         = 'data:image/jpeg;base64,' + data.image;
      resLabel.textContent = lastFilename;
      setState('result');

    } catch (err) {
      showError('Network error: ' + err.message);
      setState('idle');
    } finally {
      btnGen.disabled = false;
      btnGen.querySelector('.btn-label').textContent = 'Generate Analysis';
    }
  });

  // ── Download ──────────────────────────────────────────────────────────────
  btnDl.addEventListener('click', () => {
    if (!lastImgB64) return;
    const link    = document.createElement('a');
    link.href     = 'data:image/jpeg;base64,' + lastImgB64;
    link.download = lastFilename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  });

  // ── Allow Enter key in single inputs to trigger generate ─────────────────
  ['drug_name', 'route1_name', 'dose1', 'route2_name', 'dose2', 'mec', 'msc']
    .forEach(id => {
      const el = document.getElementById(id);
      if (el) {
        el.addEventListener('keydown', e => {
          if (e.key === 'Enter') btnGen.click();
        });
      }
    });

  // Initial state
  setState('idle');

})();
