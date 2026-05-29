@echo off
cd /d "G:\softwaropening\cyclicNoteing\webapp"
echo Starting CCS Annotation Web Server...
start "CCS-Web" python run.py
timeout /t 4 /nobreak >nul
echo Starting Cloudflare Tunnel (HTTP2)...
cloudflared.exe tunnel --url http://localhost:5000 --protocol http2
pause
