
# CNV Heatmap Generator

## Overview

This repository contains a Python script (`cnv_heatmap.py`) designed to generate a **Copy Number Variation (CNV)** heatmap from NextGene CNV tool output. The heatmap visually represents **gains**, **losses**, and **normal regions** for genes of interest, providing a clear and concise way to analyze CNV data.

The script:
1. Loads CNV data from an **Excel file**.
2. Filters and sorts the data based on **genes of interest**.
3. Saves sorted gene data to **CSV files**.
4. Creates a **heatmap PNG** to visualize gains, losses, and normal regions.

---

## Features

- **Customizable Gene Selection**: Specify any list of genes to filter and analyze.
- **Clear Visualization**: CNV gains (blue), losses (red), and normal regions (white) are easy to interpret.
- **Adaptable Paths**: File paths and sheet names are user-defined to suit various environments.
- **Flexible Output**: Save intermediate CSVs for each gene and generate a comprehensive heatmap.

---

## Requirements

- **Python Version**: 3.7+
- **Dependencies**:
  - `pandas`
  - `seaborn`
  - `matplotlib`
  - `numpy`
  - `openpyxl`

Install dependencies via pip:
```bash
pip install pandas seaborn matplotlib numpy openpyxl
```

---

## Usage

1. **Clone or Download the Repository**:
   ```bash
   git clone https://github.com/matthewlyf/cnv-heatmap-generator.git
   cd cnv-heatmap-generator
   ```

2. **Edit the Script**: Update the following variables in the `if __name__ == "__main__":` section of `cnv_heatmap.py`:
   - `EXCEL_FILE_PATH`: Path to your CNV Excel file.
   - `SHEET_NAME`: Name of the sheet in the Excel file containing CNV data.
   - `GENES_OF_INTEREST`: List of gene names you want to analyze (e.g., `["TP53", "BRCA1", "MYC"]`).
   - `OUTPUT_FOLDER`: Directory where sorted CSVs and the heatmap PNG will be saved.
   - `CSV_NAME_FOR_HEATMAP`: The specific CSV file (e.g., for a single gene) to generate the heatmap.

3. **Run the Script**:
   ```bash
   python cnv_heatmap.py
   ```

4. **Outputs**:
   - **Per-Gene CSVs**: Each gene's filtered and sorted data is saved as a separate CSV in the `OUTPUT_FOLDER`.
   - **Heatmap PNG**: The final heatmap is saved as `cnv_heatmap_with_grouped_comparisons.png` in the same folder.

---

## Example Workflow

```python
EXCEL_FILE_PATH = r"PATH/TO/YOUR/CNV_ANALYSIS.xlsx"
SHEET_NAME = "YOUR_SHEET_NAME"
GENES_OF_INTEREST = ["ATM", "BRCA1", "BRCA2", "PALB2"]
OUTPUT_FOLDER = r"PATH/TO/OUTPUT/FOLDER"
CSV_NAME_FOR_HEATMAP = "PALB2_dataframe.csv"
```

Run the script:
```bash
python cnv_heatmap.py
```

The script will:
1. Save `ATM_dataframe.csv`, `BRCA1_dataframe.csv`, etc., in the output folder.
2. Generate a heatmap PNG visualizing gains/losses.

---

## Heatmap Details

- **Color Mapping**:
  - **Blue**: CNV Gain (value >= 1.3)
  - **Red**: CNV Loss (value <= 0.7)
  - **White**: No CNV Change (0.7 < value < 1.3)

- **Features**:
  - **Column Labels**: Gene regions or descriptions (e.g., exons).
  - **Row Labels**: Sample pairs (e.g., "Sample 1" and "Sample 2").
  - **Numeric Values**: Original CNV values displayed in each heatmap cell.

---

## Folder Structure

```
cnv-heatmap-generator/
├── cnv_heatmap.py        # Main script for generating CNV heatmaps
├── README.md             # Documentation
└── requirements.txt      # Python dependencies (optional)
```

---

## Customization

- **Gene Selection**: Modify `GENES_OF_INTEREST` to analyze different genes.
- **Path Adjustments**: Update `EXCEL_FILE_PATH`, `OUTPUT_FOLDER`, and related variables to match your file system.
- **Color Mapping**: Adjust heatmap colors or thresholds by editing the `categorize` function.

---

## Notes and Limitations

1. **Input Data Format**:
   - The Excel file should include columns like `Description`, `Chr Start_x`, and sample measurement columns.
   - The script assumes a specific naming convention for sample/output columns (e.g., `"_marked_duplicates_removed_Output.pjt"`).

2. **Environment-Specific**:
   - The script includes placeholder paths (`PATH/TO/...`) that must be replaced with valid local or network paths.


---
