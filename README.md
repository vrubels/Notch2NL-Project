# Notch2NL-Project
Repository for contributing scripts and explanatory figures around Notch2NL computational workflow


This walk through starts at the conclusion of the Gordion assembler pipeline, which constructs K assembled N2NL paratypes from 10X Chromium reads. In Chimp, we see a total of 19 assembled paratypes to start, which look like this against hg38.

![assembly](figs/whole view.pdf)

The following indicates the results from this walk through, associating assembled bins with their names as they appear in our manuscript: https://doi.org/10.1101/221226

Name | Bin | Associated Genes | Genomic/BAC evidence
-----|-----|------------------|---------------------
delE2-only |	c13 | |
delE2-lrig |	c7, c16 | LRIG2 | CH251-414E24
delE2-sort	| c0, c18 | SORT1 | CH251-305M22/CH251-397N13
delE2-lrig+sort |	c10 | LRIG2, SORT1 | CH251-243D3/CH251-485H13
magi 3'	| c1, c4 | MAGI3 | AACZ03011046.1/AADA01172734.1
magi 5'	| c5, c14 | MAGI3 | AACZ03149796.1/AADA01170825.1
txnip	| c6, c8 | TXNIP | AACZ03149820.1
pde4dip	| c3, c15 | PDE4DIP | AADA01102044.1
notch2	| c2, c12, c17 | |
