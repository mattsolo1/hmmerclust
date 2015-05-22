
#hmmerclust
A python package for detection of gene clusters in bacterial genomes followed by interactive analysis of the results using ipython/jupyter notebook. It could be used for identifying genetic loci that encode biosynthetic pathways or multicomponent protein assemblies, comparing how these loci differ across taxonomical groups, and automatically extracting sequences for subsequent phylogenetic analysis, etc.

####Requirements
- python2.7/pip
- virtualenv
- hmmer installed on system
- ipython/notebook
- pandas
- matplotlib
- biopython
<hr>

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

#Installation

####Make a virtual envelope, e.g.
(see https://virtualenv.pypa.io/en/latest/)
$ brew install virtualenv
$ mkvirtualenv -p /usr/bin/python hmmerclust
$ workon hmmerclust
####Install requirements and start the notebok
$ pip install ipython[all] matplotlib pandas biopython hmmerclust$ ipython notebook
#Setting up the files
####3 key things
- Alignment folder
- Genome folder
- settings.py file


    ls

    [34malignments[m[m/             example_notebook.ipynb  [34mgenomes[m[m/                settings.py



    ls alignments/

    FlgH_PF02107_seed.txt  InvG_PF00263_seed.txt  OrgB_PB004806.txt      PrgK_PF01514_seed.txt  SpaO_PF01052_seed.txt  SpaS_PF01312_seed.txt
    InvA_PF00771_seed.txt  InvH_PF04741_seed.txt  PrgH_PF09480_seed.txt  SipB_PF04888_seed.txt  SpaP_PF00813_seed.txt
    InvC_PF00006_seed.txt  InvJ_PF02510_seed.txt  PrgI_PF09392_seed.txt  SipC_PF09599_seed.txt  SpaQ_PF01313_seed.txt
    InvE_PF07201_seed.txt  OrgA_PF09482_seed.txt  PrgJ_PB000379.txt      SipD_PF06511_seed.txt  SpaR_PF01311_seed.txt


Here are the alignment files. In this case, I queried PFAM for the proteins located in my gene cluster of interest (here, the type III secretion system). It is important to name the files starting with how you want to the search results to be annotated later, followed by an underscore. It doesn't matter what goes after the underscore.


    ls genomes/

    NC_000117.gb  NC_000907.gb  NC_000912.gb  NC_000918.gb  NC_000962.gb  NC_001264.gb  NC_002488.gb  NC_002528.gb  NC_002662.gb
    NC_000853.gb  NC_000908.gb  NC_000915.gb  NC_000919.gb  NC_000963.gb  NC_001318.gb  NC_002505.gb  NC_002570.gb
    NC_000854.gb  NC_000909.gb  NC_000916.gb  NC_000922.gb  NC_000964.gb  NC_002162.gb  NC_002506.gb  NC_002607.gb
    NC_000868.gb  NC_000911.gb  NC_000917.gb  NC_000961.gb  NC_001263.gb  NC_002163.gb  NC_002516.gb  NC_002620.gb


These are all RefSeq genbank files downloaded from NCBI. In an actual search I would use all the genomes from RefSeq. (~ 1500 chromosomes)


    # define a color dict for the proteins of interest (for locus view)
    # define an order for the heatmap 
    # if these are not defined, the will be generated automatically.
    
    !cat settings.py

    
    color_brown = 'brown'
    color_blue = 'blue'
    color_green = 'green'
    color_yellow = 'yellow'
    color_purple = 'purple'
    color_maroon = '#660033'
    color_teal = '#00FFCC'
    color_dark_grey = '#404040'
    wdomain = '#FFCC99'
    ydomain = '#FFCC66'
    chaperone = '#99CC00' #lime
    '''
    cdict = {'EccC' : color_brown,
            'EccD' : color_green,
            'MycP' : color_yellow,
            'EccB' : color_blue,
            'EccA' : color_teal,
            'EccE' : color_purple,
            'WXG100' : wdomain,
            'PE' : ydomain,
            'PPE': wdomain,
            'DUF2563' : ydomain, 
            'DUF2580' : ydomain,
            'EspA' : wdomain,
            'EspB' : wdomain,
            'EspC' : ydomain,
            'EspE' : wdomain,
            'EspF' : ydomain,
            'EspJ' : ydomain,
            'EspK' : wdomain,
            'EspL' : color_dark_grey,
            'EspD' : chaperone,
            'EspG' : chaperone,
            'EspI' : color_dark_grey,
            }
    '''
    
    hmcol = ['InvE',
            'InvC',
            'SpaP',
            'SpaQ',
            'SpaR',
            'SpaS',
            'InvA',
            'SpaO',
            'OrgB',
            'OrgA',
            'PrgK',
            'PrgH',
            'InvG',
            'InvH',
            'FlgH',
            'PrgI',
            'PrgJ',
            'InvJ',
            'SipC',
            'SipB',
            'SipD']

<hr>

#Make the database


    from hmmerclust import hmmerclust
    import settings
    
    genome_list = !ls genomes/
    genome_dir = './genomes/'
    
    db = hmmerclust.OrganismDB('t3ss_database'
                               geneome_list,
                               genome_dir,
                               freshfasta=True)

After you do this once, a fasta will be generated. If you already have the freshfasta=False, next time you make a database. The fasta can be very large; don't want to always generate this and have multiple copies on your system.

####Resultant db object structured like this:


    

<hr>

#Do the search!


    #specificy the location of the fasta 
    
    combined_fasta = './combined_fasta'
    s = hmmerclust.HmmSearch(db, combined_fasta, 
                             freshbuild=True, freshsearch=True,
                            aln_extension='.txt')

Specify the location of the combined fasta generated from the previous step and if the alignments files have an extension. Hits are added to the db object.

####Can see the organisms in the database model now have protein objects


    

#Find loci


    #args are minimum number of copies located within how many basepairs
    db.find_loci(5, 15000)


    
