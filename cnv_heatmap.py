# -*- coding: utf-8 -*-
"""
Script to generate a CNV heatmap from NextGene CNV tool output.
Paths and gene names are placeholders, so this code can be adapted
to various environments and gene lists.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


def generate_cnv_heatmap(
    excel_file_path,
    sheet_name,
    genes_of_interest,
    output_folder,
    csv_name_for_heatmap
):
    """
    Generates a CNV heatmap using data from an Excel file.

    Parameters
    ----------
    excel_file_path : str
        Path to the Excel file containing CNV data.
    sheet_name : str
        Name of the sheet in the Excel file to load.
    genes_of_interest : list of str
        A list of gene names to filter and sort.
    output_folder : str
        Path to the folder where gene-specific CSVs and the final heatmap 
        will be saved.
    csv_name_for_heatmap : str
        The specific CSV file (filtered data) used to produce the heatmap.

    Notes
    -----
    1. The script expects that the Excel data includes columns like:
       - "Description" (with gene info),
       - "Chr Start_x" (for sorting),
       - and columns indicating "output" and "sample" CNV measurements.
    2. Customize the categorization logic (`categorize`) and any columns 
       as needed to reflect your data format.
    """

    # ------------------------------------------------------------------
    # 1. LOAD EXCEL CNV DATA
    # ------------------------------------------------------------------
    df = pd.read_excel(excel_file_path, sheet_name=sheet_name)

    # Dictionary to store sorted DataFrames for each gene
    gene_dataframes = {}

    # Loop through each gene, filter rows, sort by 'Chr Start_x', then save
    for gene in genes_of_interest:
        # Filter rows where 'Description' contains the gene name
        gene_df = df[df['Description'].str.contains(gene, case=False, na=False)]

        # Example logic: reverse sorting for certain genes if needed
        # (Remove or modify if not required)
        if gene in ['BRCA1', 'PALB2']:
            sorted_gene_df = gene_df.sort_values(by='Chr Start_x', ascending=False)
        else:
            sorted_gene_df = gene_df.sort_values(by='Chr Start_x')

        # Save the sorted DataFrame to a CSV
        csv_path = f"{output_folder}/{gene}_dataframe.csv"
        sorted_gene_df.to_csv(csv_path, index=False)
        gene_dataframes[gene] = sorted_gene_df

    # ------------------------------------------------------------------
    # 2. LOAD A SPECIFIC CSV AND PREPARE HEATMAP
    # ------------------------------------------------------------------
    df_heatmap = pd.read_csv(f"{output_folder}/{csv_name_for_heatmap}")

    # Identify columns with "_marked_duplicates_removed_Output.pjt" vs. "_S"
    column_groups = {}
    for col in df_heatmap.columns:
        if "_marked_duplicates_removed_Output.pjt" in col:
            prefix = col.split("_marked_duplicates_removed_Output.pjt")[0]
            column_groups.setdefault(prefix, {}).update({"output": col})
        elif "_S" in col:
            prefix = col.split("_S")[0]
            column_groups.setdefault(prefix, {}).update({"sample": col})

    # Build pairs (output vs. sample) for each prefix
    comparison_pairs = []
    for prefix, columns in column_groups.items():
        if "output" in columns and "sample" in columns:
            comparison_pairs.append((columns["output"], columns["sample"]))

    print("Identified Comparison Pairs:", comparison_pairs)

    # If pairs are found, categorize and produce the heatmap
    if comparison_pairs:

        def categorize(value):
            """
            Categorizes numeric CNV values:
              - >= 1.3 : Gain (1)
              - <= 0.7 : Loss (-1)
              - otherwise : Normal (0)
            """
            val_rounded = round(value, 2)
            if val_rounded >= 1.3:
                return 1
            elif val_rounded <= 0.7:
                return -1
            return 0

        heatmap_data = []
        original_values = []
        row_labels = []

        # Build arrays row by row
        for (col_out, col_sample) in comparison_pairs:
            vals_out = df_heatmap[col_out]
            vals_sample = df_heatmap[col_sample]

            status_out = vals_out.apply(categorize)
            status_sample = vals_sample.apply(categorize)

            heatmap_data.append(status_out.values)
            heatmap_data.append(status_sample.values)

            original_values.append(vals_out.values)
            original_values.append(vals_sample.values)

            # Label for row axis
            row_labels.append(f"{col_out} (Sample 1)")
            row_labels.append(f"{col_sample} (Sample 2)")

        heatmap_data = np.array(heatmap_data)
        original_values = np.array(original_values)

        # Define color mapping: red (loss), white (normal), blue (gain)
        cmap = plt.cm.colors.ListedColormap(['lightcoral', 'white', 'lightblue'])
        bounds = [-1.5, -0.5, 0.5, 1.5]
        norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

        # Figure dimensions (height depends on # of row pairs)
        fig, ax = plt.subplots(figsize=(15, len(row_labels) * 0.75))  
        cax = ax.matshow(heatmap_data, cmap=cmap, norm=norm)

        # X-axis (columns) = 'Description'
        descriptions = df_heatmap['Description']
        ax.set_xticks(np.arange(len(descriptions)))
        ax.set_xticklabels(descriptions, rotation=90, ha="center", fontsize=9)

        # Y-axis (rows) = row_labels
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=9)

        # Draw cell borders + numeric values
        for (i, j), _ in np.ndenumerate(heatmap_data):
            ax.add_patch(
                patches.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="black",
                    linewidth=0.3
                )
            )
            ax.text(
                j,
                i,
                f"{original_values[i, j]:.2f}",
                ha="center",
                va="center",
                color="black",
                fontsize=7
            )

        # Thicker horizontal lines between sample pairs
        for k in range(1, len(row_labels), 2):
            ax.axhline(y=k + 0.5, color="black", linewidth=1.5)

        # Optional shading of alternate pairs for visual grouping
        for i in range(0, len(row_labels), 4):
            ax.axhspan(i - 0.5, i + 1.5, color="lightgrey", alpha=0.2)

        plt.tight_layout()

        # Save figure
        output_png = f"{output_folder}/cnv_heatmap_with_grouped_comparisons.png"
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"Heatmap saved to: {output_png}")
    else:
        print("No comparison pairs found based on column naming conventions.")


if __name__ == "__main__":

    # EXAMPLE usage: adapt these paths and lists before running
    EXCEL_FILE_PATH = r"PATH/TO/YOUR/CNV_ANALYSIS.xlsx"
    SHEET_NAME = "YOUR_SHEET"
    GENES_OF_INTEREST = ["GENE1", "GENE2", "GENE3"]  # <--- any genes of interest
    OUTPUT_FOLDER = r"PATH/TO/OUTPUT/FOLDER"
    CSV_NAME_FOR_HEATMAP = "GENE1_dataframe.csv"  # or whichever gene's CSV you want

    generate_cnv_heatmap(
        excel_file_path=EXCEL_FILE_PATH,
        sheet_name=SHEET_NAME,
        genes_of_interest=GENES_OF_INTEREST,
        output_folder=OUTPUT_FOLDER,
        csv_name_for_heatmap=CSV_NAME_FOR_HEATMAP
    )
