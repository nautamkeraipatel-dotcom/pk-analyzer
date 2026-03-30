# PK Bioavailability Analyzer 🧪

**A4-quality Pharmacokinetic Analysis Web App**  
Built by [@NautamKerai](https://github.com/NautamKerai) · GTU B.Pharm Sem 6

Generates publication-quality PK analysis charts from plasma concentration–time data.
Computes AUC, AUMC, MRT, Cmax, Tmax, Ke, t½, Cl, Vd, bioavailability (F%), loading
dose, maintenance rate, and safety status — all rendered as a 1000 DPI A4 JPEG.

---

## 📂 Project Structure

```
pk-analyzer/
├── app.py              Flask backend (routes, validation, image return)
├── pk_engine.py        PK calculation + matplotlib drawing engine
├── Bioavailability.py  Original CLI script (reference — untouched)
├── requirements.txt
├── Procfile            For Render / Railway deployment
├── templates/
│   └── index.html      Main page
└── static/
    ├── css/style.css
    └── js/main.js
```

---

## 🚀 Run Locally

```bash
# 1. Clone your repo
git clone https://github.com/<YOUR_USERNAME>/pk-analyzer.git
cd pk-analyzer

# 2. Create virtual environment
python3 -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run
python app.py

# 5. Open browser
http://localhost:5000
```

---

## ☁️ Deploy to Render (Free Hosting)

1. Push this project to a **GitHub repository**
2. Go to [render.com](https://render.com) → **New → Web Service**
3. Connect your GitHub repo
4. Set:
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `gunicorn app:app --workers 2 --timeout 120`
   - **Environment:** Python 3
5. Click **Deploy** — your site will be live at `https://your-app.onrender.com`

> **Note:** First request after inactivity may take 30–60 s (free tier spin-up).
> For production, add `gunicorn` to requirements.txt: `pip install gunicorn`

---

## 📋 How to Use

| Field | Description |
|---|---|
| Drug Name | Name of the drug (appears on chart title) |
| MEC | Minimum Effective Concentration (mg/L) — enter 0 to skip |
| MSC | Maximum Safe Concentration (mg/L) — enter 0 to skip |
| Route Name | e.g. Oral, IV, Drug A. Use "IV" in name for bioavailability formula |
| Dose (mg) | Total administered dose |
| Time Points | Space or comma-separated hours, e.g. `0 0.5 1 2 4 6 8 12 24` |
| Concentrations | Matching plasma concentrations (mg/L) |

### Example Data (Amoxicillin Oral vs IV)

**Primary Route — Oral, 500 mg**
- Time: `0 0.5 1 1.5 2 3 4 6 8 12`
- Conc: `0 1.2 3.8 5.4 5.9 4.8 3.5 1.9 0.9 0.2`

**Secondary Route — IV, 250 mg**
- Time: `0 0.25 0.5 1 2 3 4 6 8`
- Conc: `0 12.1 9.8 7.2 4.6 2.9 1.8 0.7 0.2`

---

## 🔬 PK Parameters Calculated

| Parameter | Formula |
|---|---|
| AUC | Trapezoidal rule |
| AUMC | Trapezoidal rule on C·t |
| MRT | AUMC / AUC |
| Cmax / Tmax | From spline-interpolated curve |
| Ke | \|slope\| of ln(C) vs t (terminal phase) |
| t½ | 0.693 / Ke |
| Cl | Dose / AUC |
| Vd | Cl / Ke |
| F% | (AUC₁·Dose₂) / (AUC₂·Dose₁) × 100 |
| Loading Dose | Vd × MEC |
| Maintenance Rate | Cl × MEC |

---

## 📜 License

MIT — free to use, modify, and share.
