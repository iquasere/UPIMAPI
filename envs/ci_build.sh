PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp UPIMAPI/upimapi.py UPIMAPI/resources/* "${PREFIX}/bin"
chmod +x /opt/conda/bin/upimapi.py
apt-get install -y packagekit-gtk3-module libasound2 libdbus-glib-1-2 libx11-xcb1