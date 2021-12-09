PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp UPIMAPI/upimapi.py UPIMAPI/resources/* "${PREFIX}/bin"
chmod +x /opt/conda/bin/upimapi.py