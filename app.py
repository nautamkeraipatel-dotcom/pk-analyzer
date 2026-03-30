"""
app.py — Flask backend for the PK Bioavailability Analyzer
"""

from flask import Flask, render_template, request, jsonify
import base64
import traceback
from pk_engine import generate_pk_image

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB max


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        data = request.get_json(force=True)

        # Validate required fields
        required = ['drug_name', 'mec', 'msc', 'route1_name',
                    'dose1', 'conc1_str', 'time1_str']
        for field in required:
            if not data.get(field) and data.get(field) != 0:
                return jsonify({'error': f'Missing required field: {field}'}), 400

        # Validate time/conc point counts for route 1
        try:
            t1 = [float(x) for x in data['time1_str'].replace(',', ' ').split() if x.strip()]
            c1 = [float(x) for x in data['conc1_str'].replace(',', ' ').split() if x.strip()]
        except ValueError as e:
            return jsonify({'error': f'Invalid numeric data in Route 1: {e}'}), 400

        if len(t1) < 2:
            return jsonify({'error': 'Route 1: Need at least 2 time points.'}), 400
        if len(t1) != len(c1):
            return jsonify({'error': f'Route 1: Time points ({len(t1)}) and concentrations ({len(c1)}) count must match.'}), 400

        has_second = data.get('has_second_route', False)

        if has_second:
            try:
                t2 = [float(x) for x in data.get('time2_str', '').replace(',', ' ').split() if x.strip()]
                c2 = [float(x) for x in data.get('conc2_str', '').replace(',', ' ').split() if x.strip()]
            except ValueError as e:
                return jsonify({'error': f'Invalid numeric data in Route 2: {e}'}), 400

            if len(t2) < 2:
                return jsonify({'error': 'Route 2: Need at least 2 time points.'}), 400
            if len(t2) != len(c2):
                return jsonify({'error': f'Route 2: Time points ({len(t2)}) and concentrations ({len(c2)}) count must match.'}), 400

        params = {
            'drug_name':        data['drug_name'].strip(),
            'mec':              float(data['mec']),
            'msc':              float(data['msc']),
            'route1_name':      data['route1_name'].strip(),
            'dose1':            float(data['dose1']),
            'conc1_str':        data['conc1_str'],
            'time1_str':        data['time1_str'],
            'has_second_route': has_second,
        }

        if has_second:
            params.update({
                'route2_name': data.get('route2_name', '').strip(),
                'dose2':       float(data.get('dose2', 0)),
                'conc2_str':   data.get('conc2_str', ''),
                'time2_str':   data.get('time2_str', ''),
            })

        img_bytes = generate_pk_image(params)
        img_b64   = base64.b64encode(img_bytes).decode('utf-8')

        return jsonify({
            'success': True,
            'image':   img_b64,
            'filename': f"{params['drug_name'].replace(' ', '_')}_A4_PK_Analysis.jpg"
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
