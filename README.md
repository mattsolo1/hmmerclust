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
- Python
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
First, install <a href="https://virtualenv.pypa.io/en/latest/">virtualenv</a> and <a href="https://virtualenvwrapper.readthedocs.org/en/latest/">virtualenvwrapper</a>.

Then make a virtualenv for hmmerclust:
```
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
- 16S_aligned.csv (for sorting results by phylogentic tree. Just make an empty file to start)

Before starting the search, the directory should like something like that
```
$ ls

alignments/             
example_notebook.ipynb
genomes/
settings.py
16S_aligned.csv
```
The ```alignments/``` folder should look like this:
```
$ ls alignments/

FlgH_PF02107_seed.txt  InvG_PF00263_seed.txt  OrgB_PB004806.txt      
PrgK_PF01514_seed.txt  SpaO_PF01052_seed.txt  SpaS_PF01312_seed.txt
InvA_PF00771_seed.txt  InvH_PF04741_seed.txt  PrgH_PF09480_seed.txt  
SipB_PF04888_seed.txt  SpaP_PF00813_seed.txt  InvC_PF00006_seed.txt  
InvJ_PF02510_seed.txt  PrgI_PF09392_seed.txt  SipC_PF09599_seed.txt  
SpaQ_PF01313_seed.txt  InvE_PF07201_seed.txt  OrgA_PF09482_seed.txt  
PrgJ_PB000379.txt      SipD_PF06511_seed.txt  SpaR_PF01311_seed.txt
```

In this case, I queried PFAM for the proteins located in my gene cluster of interest (here, the type III secretion system) and downloaded the hmm seed alignments. It is important to name the files starting with how you want to the search results to be annotated later, followed by an underscore. It doesn't matter what goes after the underscore.

Now, set up the genomes to query. You can use the fetch_gbwithparts() function to download genomes from NCBI.

```python
accessions_of_interest = ['NC_002516',
 'NC_002695',
 'NC_003143',
 'NC_003198',
 'NC_007645',
 'NC_007650',
 'NC_007651',
 'NC_007760',
 'NC_010287',
 'NC_013722',
 'NC_013971',
 'NC_015677',
 'NC_018518',
 'NC_022904']
 
hmmerclust.fetch_gbwithparts(accessions_of_interest, 'your@email.com', './genomes/')
 
downloading genomes... please wait
done downloading NC_002516
done downloading NC_002695
done downloading NC_003143
done downloading NC_003198
done downloading NC_007645
```

```$ !ls genomes/```
```
NC_002516.gb  NC_003198.gb  NC_007651.gb  NC_013722.gb  NC_018518.gb
NC_002695.gb  NC_007645.gb  NC_007760.gb  NC_013971.gb  NC_022904.gb
NC_003143.gb  NC_007650.gb  NC_010287.gb  NC_015677.gb
```

In an actual search I might use all the genomes from RefSeq. See a list <a href="http://www.ncbi.nlm.nih.gov/genome/browse/reference/">here</a>

In the settings file there are definitions for how results will be displayed later in some of the scripts.
The `COLOUR_DICT` variable indicates how ORFs will be coloured in the locus maps. 
The `HEATMAP_COLUMNS` variable are how proteins are ordered in the heatmap scripts. 
the `HEATMAP_ABBREVIATIONS` variable are 1 letter codes for proteins that should reflect the same order as the columns.
If these are not defined, the will be generated automatically.

```$ !cat settings.py```

```python
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

####Resultant db object
Structured so that data encapsulated in Organism objects
```python
for org in db.organisms:
    print org.name
```
```
Pseudomonas aeruginosa PAO1
Escherichia coli O157:H7 str. Sakai
Yersinia pestis CO92
Salmonella enterica subsp. enterica serovar Typhi str. CT18
Hahella chejuensis KCTC 2396
Burkholderia thailandensis E264
Burkholderia thailandensis E264
Anaeromyxobacter dehalogenans 2CP-C
Chlamydia trachomatis 434/Bu
Xanthomonas albilineans GPE PC73
Erwinia amylovora ATCC 49946
Ramlibacter tataouinensis TTB310
Bordetella pertussis 18323
Pandoraea pnomenusa 3kgm
```

<hr>

##Do the search

```python
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


####The organisms in the database model now have protein objects

```python
print db.organisms[1].proteins

NP_310807.1 - hypothetical protein,
NP_312703.1 - ATP synthase F0F1 subunit alpha,
NP_312230.1 - ABC transporter ATP-binding protein,
NP_311368.2 - glycine cleavage system transcriptional repressor,
NP_313267.1 - hypothetical protein,
NP_312611.1 - hypothetical protein,
NP_309773.1 - oligopeptide ABC transporter ATP-binding protein,
NP_309779.1 - transporter,
NP_312226.1 - hypothetical protein,
NP_311745.1 - EprI,
NP_311753.1 - surface presentation of antigens protein SpaO,
NP_312260.1 - outer membrane porin HofQ,
NP_311743.1 - EprK,
NP_311757.1 - ATP synthase SpaL,
NP_311751.1 - EpaQ,
```

##Finding loci
Now that the hits have been associated with each organism, they can be clustered into loci based on their proximity in the genome. Here, we search for loci with minimum of 5 of the search proteins within a 15k bp of each other.
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
Use the power of pandas to make sense of the results from a phylogenetic perspective, or filter the results based on the hmmer statistics.

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

Notice the `locus_id` has a Locus object, which holds information about the locus membership & location.

We can use pandas to do things like see how many loci are present in each family.
```python
# number of loci identified by family

figsize(5,5)
df.df.groupby(['org_family'])['locus_id'].nunique().plot(kind='bar')
```
<img src="/figs/locus_family.png">

...or check stats like evalue for an individual protein

```python
figsize(15,5)
df.df[df.df.hit_query=='InvG'].hit_evalue.order().plot(logy=True, kind='bar')
```

<img src="/figs/evalue.png">

##Visualize identified proteins in multiple loci and species
Plug in the pandas dataframe. If `by_locus=True`, hits shown as loci. If false, it doesn't matter where the hit is found in the genome. To order the columns in a set order, enter a list (I have this in the settings.py file). To have single letter abbreviations show up on the heatmap, set this will. To only show loci that must contain a certain protein, use the subset=[] arg. 

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
