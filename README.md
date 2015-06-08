#hmmerclust
A python package for detecting a gene cluster of interest in a set of bacterial genomes followed by interactive analysis of the results using ipython/jupyter notebook. For example, it could be used for identifying genetic loci that encode biosynthetic pathways or multicomponent protein assemblies, comparing how these loci differ across taxonomical groups, and automatically extracting sequences for subsequent phylogenetic analysis or protein expression.

Please see the <a href="https://github.com/mattsolo1/hmmerclust/tree/master/demo"><b>demo folder</b></a>, which contains <a href="https://github.com/mattsolo1/hmmerclust/blob/master/demo/hmmerclust_demo.ipynb"><b>an executable notebook demonstrating the functionality</b></a>.

####Requirements
- python2.7/pip
- virtualenv
- <a href="http://hmmer.janelia.org/software">hmmer</a> installed on system
- ipython/notebook
- pandas
- matplotlib
- biopython
<hr>

####User knowledge
- Basic python
- How to use the <a href="http://pandas.pydata.org/">pandas data analysis library</a>

####Input files
- Genbank genome sequence files 
- A folder with a set of multiple sequence alignments representing proteins of interest
- a settings file for defining colors, tables, and abbreviations
<hr>

####How it works
- Protein coding sequences are extracted from all the genome genbank files and combined into a single fasta
- Hmmbuild generates a hidden markov model from the input alignments
- Hmmsearch queries the combined proteome fasta using the hmms.
- All the data describing the organism info and taxonomy, the identified proteins, and hit statitics are stored in an object-oriented database.
- Hits are annotated based on the best-scoring hit. But other hits are preserved in the database if they need to be compared later
- The hits are clustered together based on the location of their coding sequence on the chromosome by user-specified criteria (e.g. how many CDS represent a cluster, and how far apart should they be). The hits are stored together as 'loci objects'.
- A pandas dataframe is generated with all the hit data. This allows for interactive analysis of the results, hit statistics, etc.
<hr>

####Scripts for visualizaing and analyzing results
- A heatmap for visualizing presence or absence of cluster components across a large number of species. The species can be optionally sorted based on their position in a phylogenetic tree.
- Locus maps w/ tables that describe the hits.
- Functions for extracting entire loci or sets of annotated proteins for doing multiple sequence alignments or phylogenetic analysis in a separate program.

##Installation

####Make a virtual envelope
(see https://virtualenv.pypa.io/en/latest/)
```
$ brew install virtualenv
$ mkvirtualenv -p /usr/bin/python hmmerclust
$ workon hmmerclust
```
####Install requirements and start the notebok
```
$ pip install ipython[all] matplotlib pandas biopython hmmerclust
$ ipython notebook
```
##Setting up the files
####3 key things
- Alignment folder
- Genome folder
- settings.py file

```
$ !ls

alignments/             
example_notebook.ipynb
genomes/
settings.py
```

```
$ !ls alignments/

FlgH_PF02107_seed.txt  InvG_PF00263_seed.txt  OrgB_PB004806.txt      
PrgK_PF01514_seed.txt  SpaO_PF01052_seed.txt  SpaS_PF01312_seed.txt
InvA_PF00771_seed.txt  InvH_PF04741_seed.txt  PrgH_PF09480_seed.txt  
SipB_PF04888_seed.txt  SpaP_PF00813_seed.txt  InvC_PF00006_seed.txt  
InvJ_PF02510_seed.txt  PrgI_PF09392_seed.txt  SipC_PF09599_seed.txt  
SpaQ_PF01313_seed.txt  InvE_PF07201_seed.txt  OrgA_PF09482_seed.txt  
PrgJ_PB000379.txt      SipD_PF06511_seed.txt  SpaR_PF01311_seed.txt
```

Here are the alignment files. In this case, I queried PFAM for the proteins located in my gene cluster of interest (here, the type III secretion system). It is important to name the files starting with how you want to the search results to be annotated later, followed by an underscore. It doesn't matter what goes after the underscore.

```
$ !ls genomes/

NC_000117.gb  NC_000907.gb  NC_000912.gb  NC_000918.gb  
NC_000962.gb  NC_001264.gb  NC_002488.gb  NC_002528.gb      
NC_002662.gb  NC_000853.gb  NC_000908.gb  NC_000915.gb  
NC_000919.gb  NC_000963.gb  NC_001318.gb  NC_002505.gb  
NC_002570.gb  NC_000854.gb  NC_000909.gb  NC_000916.gb  
NC_000922.gb  NC_000964.gb  NC_002162.gb  NC_002506.gb  
NC_002607.gb  NC_000868.gb  NC_000911.gb  NC_000917.gb  
NC_000961.gb  NC_001263.gb  NC_002163.gb  NC_002516.gb  
NC_002620.gb
```

These are all RefSeq genbank files downloaded from NCBI. In an actual search I would use all the genomes from RefSeq. (~ 1500 chromosomes)

In the settings file there are definitions for how results will be displayed later in some of the scripts.
The `COLOUR_DICT` variable indicates how ORFs will be coloured in the locus maps. 
The `HEATMAP_COLUMNS` variable are how proteins are ordered in the heatmap scripts. 
the `HEATMAP_ABBREVIATIONS` variable are 1 letter codes for proteins that should reflect the same order as the columns.
If these are not defined, the will be generated automatically.

$ !cat settings.py

```
COLOUR_DICT = {'InvE' : 'green', 'InvC' : 'blue'... etc}

HEATMAP_COLUMNS = ['InvE','InvC','SpaP','SpaQ','SpaS','InvA','SpaO','OrgB','OrgA','PrgK','PrgH','InvG',
        'InvH','FlgH','PrgI','PrgJ','InvJ','SipC','SipB','SipD']
        
HEATMAP_ABBREVIATIONS = ['E', 'C', 'P', 'Q', 'S', 'A', 'O', 'B', 'A', 'K', 'H', 'G',..etc]
```

<hr>

##Make the database

```python
from hmmerclust import hmmerclust
import settings

genome_list = !ls genomes/
genome_dir = './genomes/'

db = hmmerclust.OrganismDB('t3ss_database'
                           geneome_list,
                           genome_dir,
                           freshfasta=True)
```
After you do this once, a fasta will be generated. If you already have the freshfasta=False, next time you make a database. The fasta can be very large; don't want to always generate this and have multiple copies on your system.

####Resultant db object structured like this:

<hr>

##Do the search

```python
    #specificy the location of the fasta 
    
    combined_fasta = './combined_fasta'
    s = hmmerclust.HmmSearch(db, combined_fasta, 
                             freshbuild=True,
                             freshsearch=True,
                            aln_extension='.txt')
```
Specify the location of the combined fasta generated from the previous step and if the alignments files have an extension.

```
building Hmm for FlgH_PF02107_seed.txt
building Hmm for InvA_PF00771_seed.txt
building Hmm for InvC_PF00006_seed.txt
building Hmm for InvE_PF07201_seed.txt
building Hmm for InvG_PF00263_seed.txt
etc...
```
A new directory hmm/ is created. These are used in the HHsearch:

```
running HHsearch on FlgH
running HHsearch on InvA
running HHsearch on InvC
running HHsearch on InvE
running HHsearch on InvG
etc...

extracted 559 hits for FlgH.out
extracted 1000 hits for InvA.out
extracted 14722 hits for InvC.out
extracted 200 hits for InvE.out
extracted 1845 hits for InvG.out
```
Hits are added to the organisms in the db object.


####Can see the organisms in the database model now have protein objects


##Finding loci
Let's find loci with minimum of 5 of the search proteins within a 15k bp of each other
```python
db.find_loci(5, 15000)
```
```
total of 2 found for Pyrococcus abyssi GE5
finding loci for Haemophilus influenzae Rd KW20
total of 1 found for Haemophilus influenzae Rd KW20
finding loci for Mycoplasma genitalium G37
total of 1 found for Mycoplasma genitalium G37
finding loci for Methanocaldococcus jannaschii DSM 2661
total of 1 found for Methanocaldococcus jannaschii DSM 2661
finding loci for Synechocystis sp. PCC 6803
total of 0 found for Synechocystis sp. PCC 6803
etc...
```

##Make a pandas dataframe:

```python
df = hmmerclust.FinalDataFrame(db)
```

now each entry is a hit. For e.g., index the first row

```python
df.df.ix[0]
```

```
org_name                               Chlamydia trachomatis D/UW-3/CX
org_acc                                                      NC_000117
org_phylum                                                  Chlamydiae
org_class                                                 Chlamydiales
org_order                                                Chlamydiaceae
org_family                               Chlamydia/Chlamydophila group
org_genus                                                    Chlamydia
org_species                                      Chlamydia trachomatis
org_tree_order                                                       0
org_genome_length                                              1042519
org_prot_count                                                      98
org_numb_loci                                                        5
prot_acc                                                   NP_220236.1
prot_gi                                                       15605450
prot_product                    type III secretion system ATP synthase
prot_translation     MTHLQEETLLIHQWRPYRECGILSRISGSLLEAQGLSACLGELCQI...
prot_numb_of_res                                                   434
hit_query                                                         InvC
hit_evalue                                                     1.7e-51
hit_bitscore                                                     182.9
hit_bias                                                           0.1
locus_id             <hmmerclust.hmmerclust.Locus instance at 0x1af...
Name: 0, dtype: object
```

##Simultaneously visualize identified proteins in many loci at once
Plug in the pandas dataframe. If by_locus=True, hits shown as loci. If false, it doesn't matter where in the genome. To order the columns in a set order, enter a list (i have this in the settings.py file). To have single letter abbreviations show up on the heatmap, set this will. To only show loci that must contain a certain protein, use the subset=[] arg. 

```python
InvG_only = hmmerclust.HeatMap(df.df,
                  by_locus=True,
                  cols=settings.hmcol,
                  singleletters=settings.abv,
                  subset=['InvG'])
```

<img src="/figs/eg_loci.png">

##Generate a map of the locus
Any locus object can be visualized using hmmerclust.ViewLocus(locus). For example, we can generate a map of one of the loci from the InvG heatmap like so:

```python
hmmerclust.ViewLocus(InvG_only.unstacked_df.reset_index().locus_id[0])
```

<img src="/figs/eg_map.png">
