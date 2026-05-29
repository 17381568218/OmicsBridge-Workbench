"""Production launcher for CCS Annotation Web App."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')

from waitress import serve
from app import app
import socket

host_ip = socket.gethostbyname(socket.gethostname())
print(f"""
============================================
  CCS Annotation Web Server
  Local:  http://{host_ip}:5000
  (Press Ctrl+C to stop)
============================================
""")

serve(app, host='0.0.0.0', port=5000, threads=8, channel_timeout=600)
