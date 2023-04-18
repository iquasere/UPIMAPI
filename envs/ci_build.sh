PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp UPIMAPI/upimapi.py "${PREFIX}/bin"
chmod +x /opt/conda/bin/upimapi.py
ln -s "${PREFIX}/bin/upimapi.py" "${PREFIX}/bin/upimapi"