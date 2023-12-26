# MEIGN: Module Extraction for Inter-species Gene Network
This algorithm constructs a host-microbiome gene co-expression network to investigate the interactions between the host and the microbiome and identifies gene modules that closely interact. 
Specifically, it first creates three types of gene co-expression networks: host-host, host-microbiome, and microbiome-microbiome. By integrating these three networks, it constructs a host-microbiome gene co-expression network, where genes from the host and the microbiome serve as nodes, and co-expressed genes are interconnected.
Next, gene modules that contain both host and microbiome genes are identified from this gene co-expression network. This method is based on clique-search module extraction, taking into account the structure of the three networks: host-host, host-microbiome, and microbiome-microbiome. It is designed to understand interspecies interactions.

For more information on the algorithm, please see [Uehara M, Inoue T, Hase S, Sasaki E, Toyoda A, Sakakibara Y. Decoding host-microbiome interactions through co-expression network analysis within the non-human primate intestine. bioRxiv. 2023:2023-08.](https://www.biorxiv.org/content/10.1101/2023.08.11.552617v2).

## Execute the module extraction method
Execute module extraction with:
```
MEIGN.py [path of output directory] [Host-microbiome gene correlation data] [Host-host gene correlation data] [Microbiome-microbiome gene correlation data]
```

* path of output directory - set the output directory.
* Host-microbiome gene correlation data - Enter a table representing the correlation between host and microbiome genes (first column: host gene name, second column: microbiome gene name, third column: correlation coefficient). The first row will be treated as a header.
* Host-host gene correlation data - Enter a table representing the correlation between host genes (first column: host gene name, second column: host gene name, third column: correlation coefficient). The first row will be treated as a header.
* Microbiome-microbiome gene correlation data - Enter a table representing the correlation between microbiome genes (first column: microbiome gene name, second column: microbiome gene name, third column: correlation coefficient). The first row will be treated as a header.
