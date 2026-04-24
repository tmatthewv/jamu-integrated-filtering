#!/bin/bash
# ============================================================
# Batch docking: 31 ligand × 4 protein = 124 docking runs
# ============================================================

LIGAND_DIR="ligands/pdbqt"
RESULT_DIR="results"
CONFIG_DIR="config"

mkdir -p "$RESULT_DIR"

PROTEINS=("2AM9" "3ERT" "3HQR" "1SVC")
TOTAL=$(ls "$LIGAND_DIR"/*.pdbqt 2>/dev/null | wc -l)
COUNT=0
FAILED=0

echo "=============================="
echo "Batch docking — $(date)"
echo "Ligands : $TOTAL"
echo "Proteins: ${#PROTEINS[@]}"
echo "Total   : $((TOTAL * ${#PROTEINS[@]})) runs"
echo "=============================="

for PROT in "${PROTEINS[@]}"; do
    CONFIG="$CONFIG_DIR/config_${PROT}.txt"
    PROT_DIR="$RESULT_DIR/$PROT"
    mkdir -p "$PROT_DIR"

    if [ ! -f "$CONFIG" ]; then
        echo "[SKIP] Config tidak ada: $CONFIG"
        continue
    fi
    if [ ! -f "protein/${PROT}.pdbqt" ]; then
        echo "[SKIP] Protein PDBQT tidak ada: protein/${PROT}.pdbqt"
        continue
    fi

    echo ""
    echo "── Protein: $PROT ──────────────────"

    for LIGAND in "$LIGAND_DIR"/*.pdbqt; do
        NAME=$(basename "$LIGAND" .pdbqt)
        OUT="$PROT_DIR/${NAME}_out.pdbqt"
        LOG="$PROT_DIR/${NAME}_log.txt"

        # Skip jika sudah ada
        if [ -f "$LOG" ]; then
            COUNT=$((COUNT+1))
            continue
        fi

        vina \
            --config "$CONFIG" \
            --ligand  "$LIGAND" \
            --out     "$OUT" \
            --log     "$LOG" \
            2>> "$PROT_DIR/errors.txt"

        if [ $? -eq 0 ]; then
            COUNT=$((COUNT+1))
            # Ambil affinity terbaik
            AFF=$(grep -m1 "^\s*1\s" "$LOG" 2>/dev/null | awk '{print $2}')
            echo "  [${COUNT}] $NAME → $PROT : ${AFF} kcal/mol"
        else
            FAILED=$((FAILED+1))
            echo "  [FAIL] $NAME → $PROT"
        fi
    done
done

echo ""
echo "=============================="
echo "Selesai: $COUNT berhasil, $FAILED gagal"
echo "Waktu  : $(date)"
echo "=============================="
