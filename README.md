# HisTrader
A Tool to Identify Nucleosome Free Regions from ChIP-Seq Against Histone Modifications

USAGE: perl HISTRADER.pl --bedGraph ChIP.bedGraph --peaks ChIP.bed

OPTIONS:

--bedGraph    Specify ChIP-Seq signal file (bedGraph format) (required)

--peaks		    Specify the BROAD PEAK file (bed format) (required)

--genome	    Specify genome fasta file (optional)
		          Used to extract the DNA sequence within the valleys/footprints/NFRs

--trim		    Specify that fasta sequences should be trimmed
		          Default (OFF).	REQUIRES --genome and --trimSize.

--trimSize	  The length of fasta sequence in base pairs (bp) returned for each valley / footprint / NFR. Sequences are     trimmed from the centre of each valley / footprint / NFR.
		          REQUIRES --trim and --genome

--out		      The output file prefix
		          Default=Histrader (optional)

--method    	Method used to call valleys / footprints/ NFRs
	          	Options = MA (moving average) / DIFF (differencing) / BOTH (merge of MA and DIFF)
		          Default = BOTH

--step	    	The fixed step size to use. The ChIP-seq profile will be converted to have fixed step equal to this number in   base pairs (bp).
	          	Default=25 (optional)

--minSize	    The mininum peak size. Valleys/Footprints/NFRs will NOT be called in peaks that are smaller than this value
		          Default=500 (optional)

--nucSize   	The estimated nucleosome size in base pairs (bp).This value is used for the moving average and smoothing the signal.
	           	This value should be divisible by the step size specified (--step).
	          	Default: 150 (~147).	150 / 25 = 6 (6 bins is equivalent to 1 nucleosome)

--mergeMulti	The multiplier of step (--step) to use for merging.
		          Default = 3	( 3 x 25 = 75 or ~ 1/2 a nucleosome)

--maMulti	    The moving average multiplier for the slow moving average.
		          Default = 3 (ie. Slow moving average is equilavent to 3 nucleosomes [3 x 150])

--outBG	    	Output the fixed step ChIP-Seq signal at the specified broad peaks (bedGraph format).

--help        (prints this message)
