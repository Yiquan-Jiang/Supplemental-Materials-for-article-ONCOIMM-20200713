#  Investigating Mechanisms of Response or Resistance to Immune Checkpoint Inhibitors by Analyzing Cell-Cell Communications in Tumors Before and After Programmed Cell Death-1 (PD-1) Targeted Therapy: An Integrative Analysis Using Single-cell RNA and Bulk-RNA Sequencing Data

Background: Currently, a significant proportion of cancer patients do not benefit from programmed cell death-1 (PD-1)-targeted therapy. Overcoming drug resistance remains a challenge.

Methods: Single-cell RNA sequencing data from 53,030 cells collected before and after anti-PD1 therapy were analyzed. Cell-cell interaction analyses were performed to determine the differences between pretreatment responders and nonresponders and the relative differences in changes from pretreatment to posttreatment status between responders and nonresponders to ultimately investigate the specific mechanisms underlying response and resistance to anti-PD1 therapy. Bulk-RNA sequencing data were used to validate our results. Furthermore, we analyzed the evolutionary trajectory of ligands/receptors in specific cell types in responders and nonresponders.

Results: We identified several different cell-cell interactions, like WNT5A-PTPRK, EGFR-AREG, AXL-GAS6 and ACKR3-CXCL12, based on pretreatment data from responders and nonresponders. Furthermore, relative differences in the changes from pretreatment to posttreatment status between responders and nonresponders existed in SELE-PSGL-1, CXCR3-CCL19, CCL4-SLC7A1, CXCL12-CXCR3, EGFR-AREG, THBS1-a3b1 complex, TNF-TNFRSF1A, TNF-FAS and TNFSF10-TNFRSF10D interactions. In trajectory analyses of tumor-specific exhausted CD8 T cells using ligand/receptor genes, we identified a cluster of T cells that presented a distinct pattern of ligand/receptor expression. They highly expressed suppressive genes like HAVCR2 and KLRC1, cytotoxic genes like GZMB and FASLG and the tissue-residence-related gene CCL5. These cells had increased expression of survival-related and tissue-residence-related genes, like heat shock protein genes and the interleukin-7 receptor (IL-7R), CACYBP and IFITM3 genes, after anti-PD1 therapy.

Conclusions: These results reveal the mechanisms underlying anti-PD1 therapy response and offer abundant clues for potential strategies to improve immunotherapy.


# Install CellPhoneDB V2
Please download and install the CellPhoneDB V2 at https://github.com/Teichlab/cellphonedb.



# Supplemental-Materials-for-article-ONCOIMM-20200713
The r codes used for data analyses and figure formation are deposited in The_core_R_Script.R

Supplemental Materials list:

Supplemental Table 1 Sample Information

Supplemental Table 2 Ligand-Receptor Pair List and Ligand Receptor Genes

Supplemental Table 3 Different Cell-Cell Interactions Between Pretreatment Responders and Nonresponders

Supplemental Table 4 Relative Differences in Changes from Pretreatment to Posttreatment Status Between Responders and Nonresponders

Supplemental Table 5 Differentially Expressed Genes Between CAFs and Myofibroblasts	

Supplemental Table 6: GSVA between CAFs and Myofibroblasts

Supplemental Table 7 Differentially Expressed Ligands/Receptors Across the 4 Branches of Exhausted CD8 T cells in Responders

Supplemental Table 8 GSVA of the 4 Branches of Exhausted CD8 T cells in Responders

Supplemental Table 9 Differentially Expressed Genes Between Pretreatment- and Posttreatment Status in Cluster 1 Exhausted CD8 T cells in responders

Fig. S1 Specific Comparison of Each Ligand-receptor Interaction in Pretreatment responders with Each Ligand-receptor Interaction in Pretreatment Nonresponders. Notes, a Ratio (pretreatment responders/nonresponders) >1 indicates that a higher interaction intensity existed in responders than in nonresponders (shown in red). A Ratio (pretreatment responders/nonresponders) <1 indicates that a lower interaction intensity exists in responders than in nonresponders (shown in blue).

Fig. S2 Overlapping Genes Between Our Ligand or Receptor Genes and the 693 Differentially Expressed Genes (DEGs) Identified in the Hugo et al. Study Between Pretreatment Responders and Pretreatment Nonresponders. Notes, a Ratio (pretreatment responders/nonresponders) >1 indicates that a higher interaction intensity existed in responders than in nonresponders (shown in red). A Ratio (pretreatment responders/nonresponders) <1 indicates that a lower interaction intensity exists in responders than in nonresponders (shown in blue).

Fig. S3 Relative Differences (Responders vs. Nonresponders) in Changes from Pretreatment to Posttreatment Status. Notes, a “Relative Ratio” between responders and nonresponders >1 means that the interaction intensity was relatively increased in responders or relatively decreased in nonresponders during treatment (shown in red). A “relative ratio” between responders and nonresponders <1 means that the interaction intensity was relatively decreased in responders or relatively increased in nonresponders during treatment (shown in blue).

Fig.S4 Overlapping Genes Between Our Ligand or Receptor Genes and with 2,670 DEGs that changed differentially from pretreatment to posttreatment status between responders and nonresponders in the Riaz et al. study. Notes, a “Relative Ratio” between responders and nonresponders >1 means that the interaction intensity was relatively increased in responders or relatively decreased in nonresponders during treatment (shown in red). A “relative ratio” between responders and nonresponders <1 means that the interaction intensity was relatively decreased in responders or relatively increased in nonresponders during treatment (shown in blue).

Fig. S5 DEGs and GSVA between CAFs and Myofibroblasts.

Fig. S6 Validation of the Ligand/receptor Expression Patterns Identified in Our Study with Other Immunotherapy Datasets using the TIDE Platform. Notes, genes that showed all red or all blue indicated a high correlation with immunotherapy response.
