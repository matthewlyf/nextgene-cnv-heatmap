
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


"""
Main function to load CNV data from NextGene output and produce a heatmap
indicating gains/losses by gene or region.
"""
# ----------------------------------------
# 1. Load the main Excel file
# ----------------------------------------
file_path = r"Y:/DataAnalysis/Genome Analysts/ML/CNV analysis/CNV_ANALYSIS_DUPLICATE_READ_REMOVAL_28OCT2024.xlsx"
sheet_name = "22SOMATICNGS6_NO_CONTROLS"
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Genes of interest
genes_of_interest = ['ATM', 'PALB2', 'BRCA1', 'BRCA2']
gene_dataframes = {}

# Sort the data by each gene (reverse for BRCA1/PALB2, normal for others)
for gene in genes_of_interest:
    gene_df = df[df['Description'].str.contains(gene, case=False, na=False)]
    if gene in ['BRCA1', 'PALB2']:
        sorted_gene_df = gene_df.sort_values(by='Chr Start_x', ascending=False)
    else:
        sorted_gene_df = gene_df.sort_values(by='Chr Start_x')

    # Save as CSV (optional)
    sorted_gene_df.to_csv(f"Y:/DataAnalysis/Genome Analysts/ML/CNV analysis/{gene}_dataframe.csv", index=False)
    gene_dataframes[gene] = sorted_gene_df

# ----------------------------------------
# 2. Load a specific CSV for generating heatmap
# ----------------------------------------
df_heatmap = pd.read_csv(r"Y:\DataAnalysis\Genome Analysts\ML\CNV analysis\22SOMATICNGS6_No_Controls\PALB2_dataframe.csv")

# Identify comparison pairs (col1: Output, col2: Sample)
column_groups = {}
for col in df_heatmap.columns:
    if "_marked_duplicates_removed_Output.pjt" in col:
        prefix = col.split("_marked_duplicates_removed_Output.pjt")[0]
        column_groups.setdefault(prefix, {}).update({"output": col})
    elif "_S" in col:
        prefix = col.split("_S")[0]
        column_groups.setdefault(prefix, {}).update({"sample": col})

comparison_pairs = []
for prefix, columns in column_groups.items():
    if "output" in columns and "sample" in columns:
        comparison_pairs.append((columns["output"], columns["sample"]))

# Debug / check
print("Identified Comparison Pairs:", comparison_pairs)

# ----------------------------------------
# 3. Create and save the heatmap
# ----------------------------------------
if comparison_pairs:
    def categorize(value):
        val_rounded = round(value, 2)
        if val_rounded >= 1.3:
            return 1  # Gain
        elif val_rounded <= 0.7:
            return -1  # Loss
        return 0      # Normal

    heatmap_data = []
    original_values = []
    row_labels = []

    for (col1, col2) in comparison_pairs:
        values1 = df_heatmap[col1]
        values2 = df_heatmap[col2]

        status1 = values1.apply(categorize)
        status2 = values2.apply(categorize)

        heatmap_data.append(status1.values)
        heatmap_data.append(status2.values)

        original_values.append(values1.values)
        original_values.append(values2.values)

        row_labels.append(f"{col1} (Sample 1)")
        row_labels.append(f"{col2} (Sample 2)")

    heatmap_data = np.array(heatmap_data)
    original_values = np.array(original_values)

    # Light Red (-1), White (0), Light Blue (1)
    cmap = plt.cm.colors.ListedColormap(['lightcoral', 'white', 'lightblue'])
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(15, len(comparison_pairs) * 1.5))
    cax = ax.matshow(heatmap_data, cmap=cmap, norm=norm)

    # X-axis: Exon Descriptions
    ax.set_xticks(np.arange(len(df_heatmap['Description'])))
    ax.set_xticklabels(df_heatmap['Description'], rotation=90, ha="center", fontsize=10)

    # Y-axis: Sample Labels
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)

    # Add black borders + numeric values
    for (i, j), val in np.ndenumerate(heatmap_data):
        ax.add_patch(
            patches.Rectangle(
                (j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", linewidth=0.5
            )
        )
        ax.text(j, i, f"{original_values[i, j]:.2f}", ha="center", va="center", color="black", fontsize=8)

    # Thicker line between each sample pair
    for k in range(1, len(row_labels), 2):
        ax.axhline(y=k + 0.5, color="black", linewidth=3)

    # Shading for alternate sample pairs
    for i in range(0, len(row_labels), 4):
        ax.axhspan(i - 0.5, i + 1.5, color="lightgrey", alpha=0.2)

    plt.tight_layout()
    output_png = r"Y:\DataAnalysis\Genome Analysts\ML\CNV analysis\22SOMATICNGS6_No_Controls\cnv_heatmap_with_grouped_comparisons_no_legend.png"
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Heatmap saved to: {output_png}")
else:
    print("No comparison pairs found based on the naming convention.")

