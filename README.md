{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#hmmerclust\n",
    "A python package for detection of gene clusters in bacterial genomes followed by interactive analysis of the results using ipython/jupyter notebook. It could be used for identifying genetic loci that encode biosynthetic pathways or multicomponent protein assemblies, comparing how these loci differ across taxonomical groups, and automatically extracting sequences for subsequent phylogenetic analysis, etc.\n",
    "\n",
    "####Requirements\n",
    "- python2.7/pip\n",
    "- virtualenv\n",
    "- hmmer installed on system\n",
    "- ipython/notebook\n",
    "- pandas\n",
    "- matplotlib\n",
    "- biopython\n",
    "<hr>\n",
    "\n",
    "####Input files\n",
    "- Genbank genome sequence files \n",
    "- A folder with a set of multiple sequence alignments representing proteins of interest\n",
    "- a settings file for defining colors, tables, and abbreviations\n",
    "<hr>\n",
    "\n",
    "####How it works\n",
    "- Protein coding sequences are extracted from all the genome genbank files and combined into a single fasta\n",
    "- Hmmbuild generates a hidden markov model from the input alignments\n",
    "- Hmmsearch queries the combined proteome fasta using the hmms.\n",
    "- All the data describing the organism info and taxonomy, the identified proteins, and hit statitics are stored in an object-oriented database.\n",
    "- Hits are annotated based on the best-scoring hit. But other hits are preserved in the database if they need to be compared later\n",
    "- The hits are clustered together based on the location of their coding sequence on the chromosome by user-specified criteria (e.g. how many CDS represent a cluster, and how far apart should they be). The hits are stored together as 'loci objects'.\n",
    "- A pandas dataframe is generated with all the hit data. This allows for interactive analysis of the results, hit statistics, etc.\n",
    "<hr>\n",
    "\n",
    "####Scripts for visualizaing and analyzing results\n",
    "- A heatmap for visualizing presence or absence of cluster components across a large number of species. The species can be optionally sorted based on their position in a phylogenetic tree.\n",
    "- Locus maps w/ tables that describe the hits.\n",
    "- Functions for extracting entire loci or sets of annotated proteins for doing multiple sequence alignments or phylogenetic analysis in a separate program.\n",
    "\n",
    "#Installation\n",
    "\n",
    "####Make a virtual envelope, e.g.\n",
    "(see https://virtualenv.pypa.io/en/latest/)"
   ]
  },
  {
   "cell_type": "raw",
   "metadata": {},
   "source": [
    "$ brew install virtualenv\n",
    "$ mkvirtualenv -p /usr/bin/python hmmerclust\n",
    "$ workon hmmerclust"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "####Install requirements and start the notebok"
   ]
  },
  {
   "cell_type": "raw",
   "metadata": {},
   "source": [
    "$ pip install ipython[all] matplotlib pandas biopython hmmerclust"
   ]
  },
  {
   "cell_type": "raw",
   "metadata": {},
   "source": [
    "$ ipython notebook"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#Setting up the files\n",
    "####3 key things\n",
    "- Alignment folder\n",
    "- Genome folder\n",
    "- settings.py file"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {
    "collapsed": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\u001b[34malignments\u001b[m\u001b[m/             example_notebook.ipynb  \u001b[34mgenomes\u001b[m\u001b[m/                settings.py\r\n"
     ]
    }
   ],
   "source": [
    "ls"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {
    "collapsed": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "FlgH_PF02107_seed.txt  InvG_PF00263_seed.txt  OrgB_PB004806.txt      PrgK_PF01514_seed.txt  SpaO_PF01052_seed.txt  SpaS_PF01312_seed.txt\r\n",
      "InvA_PF00771_seed.txt  InvH_PF04741_seed.txt  PrgH_PF09480_seed.txt  SipB_PF04888_seed.txt  SpaP_PF00813_seed.txt\r\n",
      "InvC_PF00006_seed.txt  InvJ_PF02510_seed.txt  PrgI_PF09392_seed.txt  SipC_PF09599_seed.txt  SpaQ_PF01313_seed.txt\r\n",
      "InvE_PF07201_seed.txt  OrgA_PF09482_seed.txt  PrgJ_PB000379.txt      SipD_PF06511_seed.txt  SpaR_PF01311_seed.txt\r\n"
     ]
    }
   ],
   "source": [
    "ls alignments/"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Here are the alignment files. In this case, I queried PFAM for the proteins located in my gene cluster of interest (here, the type III secretion system). It is important to name the files starting with how you want to the search results to be annotated later, followed by an underscore. It doesn't matter what goes after the underscore."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {
    "collapsed": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "NC_000117.gb  NC_000907.gb  NC_000912.gb  NC_000918.gb  NC_000962.gb  NC_001264.gb  NC_002488.gb  NC_002528.gb  NC_002662.gb\r\n",
      "NC_000853.gb  NC_000908.gb  NC_000915.gb  NC_000919.gb  NC_000963.gb  NC_001318.gb  NC_002505.gb  NC_002570.gb\r\n",
      "NC_000854.gb  NC_000909.gb  NC_000916.gb  NC_000922.gb  NC_000964.gb  NC_002162.gb  NC_002506.gb  NC_002607.gb\r\n",
      "NC_000868.gb  NC_000911.gb  NC_000917.gb  NC_000961.gb  NC_001263.gb  NC_002163.gb  NC_002516.gb  NC_002620.gb\r\n"
     ]
    }
   ],
   "source": [
    "ls genomes/"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "These are all RefSeq genbank files downloaded from NCBI. In an actual search I would use all the genomes from RefSeq. (~ 1500 chromosomes)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {
    "collapsed": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\r\n",
      "color_brown = 'brown'\r\n",
      "color_blue = 'blue'\r\n",
      "color_green = 'green'\r\n",
      "color_yellow = 'yellow'\r\n",
      "color_purple = 'purple'\r\n",
      "color_maroon = '#660033'\r\n",
      "color_teal = '#00FFCC'\r\n",
      "color_dark_grey = '#404040'\r\n",
      "wdomain = '#FFCC99'\r\n",
      "ydomain = '#FFCC66'\r\n",
      "chaperone = '#99CC00' #lime\r\n",
      "'''\r\n",
      "cdict = {'EccC' : color_brown,\r\n",
      "        'EccD' : color_green,\r\n",
      "        'MycP' : color_yellow,\r\n",
      "        'EccB' : color_blue,\r\n",
      "        'EccA' : color_teal,\r\n",
      "        'EccE' : color_purple,\r\n",
      "        'WXG100' : wdomain,\r\n",
      "        'PE' : ydomain,\r\n",
      "        'PPE': wdomain,\r\n",
      "        'DUF2563' : ydomain, \r\n",
      "        'DUF2580' : ydomain,\r\n",
      "        'EspA' : wdomain,\r\n",
      "        'EspB' : wdomain,\r\n",
      "        'EspC' : ydomain,\r\n",
      "        'EspE' : wdomain,\r\n",
      "        'EspF' : ydomain,\r\n",
      "        'EspJ' : ydomain,\r\n",
      "        'EspK' : wdomain,\r\n",
      "        'EspL' : color_dark_grey,\r\n",
      "        'EspD' : chaperone,\r\n",
      "        'EspG' : chaperone,\r\n",
      "        'EspI' : color_dark_grey,\r\n",
      "        }\r\n",
      "'''\r\n",
      "\r\n",
      "hmcol = ['InvE',\r\n",
      "        'InvC',\r\n",
      "        'SpaP',\r\n",
      "        'SpaQ',\r\n",
      "        'SpaR',\r\n",
      "        'SpaS',\r\n",
      "        'InvA',\r\n",
      "        'SpaO',\r\n",
      "        'OrgB',\r\n",
      "        'OrgA',\r\n",
      "        'PrgK',\r\n",
      "        'PrgH',\r\n",
      "        'InvG',\r\n",
      "        'InvH',\r\n",
      "        'FlgH',\r\n",
      "        'PrgI',\r\n",
      "        'PrgJ',\r\n",
      "        'InvJ',\r\n",
      "        'SipC',\r\n",
      "        'SipB',\r\n",
      "        'SipD']"
     ]
    }
   ],
   "source": [
    "# define a color dict for the proteins of interest (for locus view)\n",
    "# define an order for the heatmap \n",
    "# if these are not defined, the will be generated automatically.\n",
    "\n",
    "!cat settings.py"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "<hr>"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#Make the database"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "from hmmerclust import hmmerclust\n",
    "import settings\n",
    "\n",
    "genome_list = !ls genomes/\n",
    "genome_dir = './genomes/'\n",
    "\n",
    "db = hmmerclust.OrganismDB('t3ss_database'\n",
    "                           geneome_list,\n",
    "                           genome_dir,\n",
    "                           freshfasta=True)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "After you do this once, a fasta will be generated. If you already have the freshfasta=False, next time you make a database. The fasta can be very large; don't want to always generate this and have multiple copies on your system."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "####Resultant db object structured like this:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "<hr>"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#Do the search!"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "#specificy the location of the fasta \n",
    "\n",
    "combined_fasta = './combined_fasta'\n",
    "s = hmmerclust.HmmSearch(db, combined_fasta, \n",
    "                         freshbuild=True, freshsearch=True,\n",
    "                        aln_extension='.txt')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Specify the location of the combined fasta generated from the previous step and if the alignments files have an extension. Hits are added to the db object."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "####Can see the organisms in the database model now have protein objects"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#Find loci"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "#args are minimum number of copies located within how many basepairs\n",
    "db.find_loci(5, 15000)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 2",
   "language": "python",
   "name": "python2"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 0
}
