PREFIX="/opt/conda"
mkdir -p "${PREFIX}/bin"
cp UPIMAPI/upimapi.py UPIMAPI/uniprot_support.py "${PREFIX}/bin"
chmod +x /opt/conda/bin/upimapi.py