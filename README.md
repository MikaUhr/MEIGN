# MEIGN: Module Extraction for Inter-species Gene Network
This algorithm constructs a host-microbiome gene co-expression network to investigate the interactions between the host and the microbiome and identifies gene modules that closely interact. 
Specifically, it first creates three types of gene co-expression networks: host-host, host-microbiome, and microbiome-microbiome. By integrating these three networks, it constructs a host-microbiome gene co-expression network, where genes from the host and the microbiome serve as nodes, and co-expressed genes are interconnected.
Next, gene modules that contain both host and microbiome genes are identified from this gene co-expression network. This method is based on clique-search module extraction, taking into account the structure of the three networks: host-host, host-microbiome, and microbiome-microbiome. It is designed to understand interspecies interactions.

For more information on the algorithm, please see [Uehara M, Inoue T, Hase S, Sasaki E, Toyoda A, Sakakibara Y. Decoding host-microbiome interactions through co-expression network analysis within the non-human primate intestine. bioRxiv. 2023:2023-08.](https://www.biorxiv.org/content/10.1101/2023.08.11.552617v2).


## Execute the module extraction method
Execute module extraction with:
```
MEIGN.py [Path of output directory] [Host-microbiome gene correlation data] [Host-host gene correlation data] [Microbiome-microbiome gene correlation data]
```

* Path of output directory - set the output directory.
* Host-microbiome gene correlation data - Enter a tab-delimited table representing the correlation between host and microbiome genes (first column: host gene name, second column: microbiome gene name, third column: correlation coefficient). The first row will be treated as a header.
* Host-host gene correlation data - Enter a tab-delimited table representing the correlation between host genes (first column: host gene name, second column: host gene name, third column: correlation coefficient). The first row will be treated as a header.
* Microbiome-microbiome gene correlation data - Enter a tab-delimited table representing the correlation between microbiome genes (first column: microbiome gene name, second column: microbiome gene name, third column: correlation coefficient). The first row will be treated as a header.

### Option
```
MEIGN.py [Path of output directory] [Host-microbiome gene correlation data] [Host-host gene correlation data] [Microbiome-microbiome gene correlation data] --edge [float value] --number [int value] --cliqueness [float value]
```
* --edge, -e [float value] : Correlation coefficient threshold. An edge is drawn if this threshold is exceeded in the gene co-expression network. 0.0 - 1.0 is accepted.
* --number, -n [int value]: Minimum number of nodes in clique.
* --cliqueness, -c [float value]: Threshold of cliqueness when merging modules. 0.0 - 1.0 is accepted.


## Example
### Input files

Host-microbiome gene correlation data
```
Host_Gene  Microbiome_Gene  Correlation_Coefficient
Gene_h1  Gene_m1  0.480
Gene_h1  Gene_m2  -0.611
Gene_h1  Gene_m3  0.252
Gene_h1  Gene_m4  0.552
Gene_h2  Gene_m1  -0.840
Gene_h2  Gene_m2  0.601
Gene_h2  Gene_m3  0.194
Gene_h2  Gene_m4  0.941
Gene_h3  Gene_m1  0.940
Gene_h3  Gene_m2  0.765
Gene_h3  Gene_m3  -0.522
Gene_h3  Gene_m4  0.221
Gene_h4  Gene_m1  0.940
Gene_h4  Gene_m2  0.765
Gene_h4  Gene_m3  -0.522
Gene_h4  Gene_m4  -0.429
```
Host-host gene correlation data
```
Host_Gene_1  Host_Gene_2  Correlation_Coefficient
Gene_h1  Gene_h2  0.840
Gene_h1  Gene_h3  0.601
Gene_h2  Gene_h3  -0.957
Gene_h2 Gene_h4  0.943
Gene_h3 Gene_h4  -0.431
```

Microbiome-microbiome gene correlation data
```
Microbiome_Gene_1  Microbiome_Gene_2  Correlation_Coefficient
Gene_m1  Gene_m2  0.134
Gene_m1  Gene_m3  -0.721
Gene_m1  Gene_m4  0.922
Gene_m2  Gene_m3  0.957
Gene_m2 Gene_m4  -0.434
Gene_m3 Gene_m4  0.831
```

### Output file
A file named "module_list.tsv" is output to the specified output directory. This file represents a single module on a single line, and the genes contained in the module are described in a tab-delimited list.
```
Gene_h1  Gene_h2  Gene_h3  Gene_m1
Gene_h4  Gene_m2  Gene_m3
```
