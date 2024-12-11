.PHONY: freeze
freeze:
	conda env export | grep -v prefix > conda.env.yml

.PHONY: osx-arm
osx-arm:
	CONDA_SUBDIR=osx-64 conda create -n rnaseqdiff
	conda activate rnaseqdiff
	conda config --env --set subdir osx-64
	conda config --set channel_priority flexible

.PHONY: PerlBio
PerlBio:
	cpanm -n Bio::Perl
