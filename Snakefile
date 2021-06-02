configfile: "config.json"

rule all:
	input:
		expand("hyphy/{locus}/done.txt", locus=config["locus"])
	
rule translate:
	input:
		"genes/{locus}/"
	output:
		"genes/{locus}/codon_guided_msa/translated.fasta"
	params:
		"genes/{locus}/*fasta"
	threads: 1

	shell:
		"ls {params} | xargs -I@ -n 1 sh -c 'seqkit translate -M --clean @ >> {output}'"
configfile: "config.json"

rule remove_stop:
	input:
		"genes/{locus}/codon_guided_msa/translated.fasta"
	output:
		"genes/{locus}/codon_guided_msa/translated_nostop.fasta"
	shell:
		"scripts/remove_stop.py {input} {output}"


rule protein_msa:
	input:
		"genes/{locus}/codon_guided_msa/translated_nostop.fasta"
	output:
		"genes/{locus}/codon_guided_msa/msa_translated.fasta"
	params:
		"--maxiterate 1000 --localpair"
	conda:
		"envs/pos_selection.yaml"

	threads: 1
	shell:
		"mafft --thread 1 --quiet {params} {input} > {output}"

	
rule merge_cds:
	input:
		"genes/{locus}"
	output:
		"genes/{locus}/codon_guided_msa/multi_fasta_cds.fasta"
	params:
		"genes/{locus}/*fasta"
	threads: 1

	shell:
		"ls {params} | xargs -I@ sh -c 'cat @ >> {output}'"

rule codon_guided_msa:
	input:
		cds_fasta="genes/{locus}/codon_guided_msa/multi_fasta_cds.fasta",
		prot_alnm="genes/{locus}/codon_guided_msa/msa_translated.fasta"
	output:
		"genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"
	conda:
		"envs/pos_selection.yaml"
	threads: 1

	shell:
		"pal2nal.pl {input.prot_alnm} {input.cds_fasta} -output fasta > {output}"

rule raxml:
	input:
		"genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"
	output:
		besttree="RAxML_bestTree.{locus}_raxml",
		info="RAxML_info.{locus}_raxml",
		log_file="RAxML_log.{locus}_raxml",
		pars="RAxML_parsimonyTree.{locus}_raxml",
		result="RAxML_result.{locus}_raxml"
	params:
		model="GTRGAMMA",
		suffix="{locus}_raxml"

	threads: 1
	conda:
		"envs/pos_selection.yaml"
	shell:
		"raxmlHPC -T 1 -s {input} -n {params.suffix} -m {params.model} -p 123"
		

rule tidy_up:
	input:
		besttree="RAxML_bestTree.{locus}_raxml",
		info="RAxML_info.{locus}_raxml",
		log_file="RAxML_log.{locus}_raxml",
		pars="RAxML_parsimonyTree.{locus}_raxml",
		result="RAxML_result.{locus}_raxml"
	output:
		besttree="raxml/{locus}/RAxML_bestTree.{locus}_raxml",
		info="raxml/{locus}/RAxML_info.{locus}_raxml",
		log_file="raxml/{locus}/RAxML_log.{locus}_raxml",
		pars="raxml/{locus}/RAxML_parsimonyTree.{locus}_raxml",
		result="raxml/{locus}/RAxML_result.{locus}_raxml"

	threads: 1
		
	shell:
		"mv {input.besttree} {output.besttree}; mv {input.info} {output.info}; mv {input.log_file} {output.log_file}; mv {input.pars} {output.pars}; mv {input.result} {output.result}"

rule hyphy_fubar:
	input:
		besttree="raxml/{locus}/RAxML_bestTree.{locus}_raxml",
		msa="genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"
	output:
		summary="hyphy/{locus}/fubar/{locus}.FUBAR.summary"
	conda:
		"envs/hyphy.yaml"

	threads: 1

	shell:
		"hyphy fubar --alignment {input.msa} --tree {input.besttree} > {output.summary}"


rule hyphy_fel:
	input:
		besttree="raxml/{locus}/RAxML_bestTree.{locus}_raxml",
		msa="genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"
	output:
		summary="hyphy/{locus}/fel/{locus}.FEL.summary"
	conda:
		"envs/hyphy.yaml"

	threads: 1

	shell:
		"hyphy fel --alignment {input.msa} --tree {input.besttree} > {output.summary}"

rule codeml_setup:
	input:
		besttree="raxml/{locus}/RAxML_bestTree.{locus}_raxml",
		msa="genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"
	output:
		directory=directory("codeml/{locus}/ctl_files")
		
	params:
		locus="{locus}",
		besttree="../../raxml/{locus}/RAxML_bestTree.{locus}_raxml",
		msa="../../genes/{locus}/codon_guided_msa/{locus}_codon_guided_msa.fasta"

	conda:
		"envs/pos_selection.yaml"
	threads: 1

	shell:
		"mkdir {output.directory};"
		"scripts/setup_codeml_ctl_file.py {params.msa} {params.besttree} {params.locus}"

rule codeml_run:
	input:
		directory="codeml/{locus}/ctl_files"
	output:
		directory=directory("codeml/{locus}/mlc_files"),
		logs=directory("codeml/{locus}/log_files")
	params:
		wdir="codeml/{locus}/",

		one="ctl_files/{locus}_1_codeml.ctl",
		one_log="log_files/{locus}_1_codeml.log",

		two="ctl_files/{locus}_2_codeml.ctl",
		two_log="log_files/{locus}_2_codeml.log",

		seven="ctl_files/{locus}_7_codeml.ctl",
		seven_log="log_files/{locus}_7_codeml.log",

		eight="ctl_files/{locus}_8_codeml.ctl",
		eight_log="log_files/{locus}_8_codeml.log",

		tidy = "rm 2NG*; rm lnf; rm rst*; rm rub"


	threads: 1
	conda:
		"envs/pos_selection.yaml"
	shell:
		"mkdir {output.directory};"
		"mkdir {output.logs};"
		"cd {params.wdir};"
		"codeml {params.one} > {params.one_log};"
		"codeml {params.two} > {params.two_log};"
		"codeml {params.seven} > {params.seven_log};"
		"codeml {params.eight} > {params.eight_log};"
		"{params.tidy}"

rule codeml_stats:
	input:
		"codeml/{locus}/mlc_files/"
	output:
		directory=directory("codeml/{locus}/stats")
	params:
		mlcs="codeml/{locus}/mlc_files/mlc*",
		stats="codeml/{locus}/stats/{locus}_modeltesting.txt",
		summary_models="codeml/{locus}/stats/{locus}_summary_models.txt",
		summary_stats="codeml/{locus}/stats/{locus}_summary_stats.txt"

	threads: 1

	shell:
		"mkdir {output.directory};"
		"grep lnL {params.mlcs} > {params.stats};"
		"scripts/codeml_stats.py {params.stats} {params.summary_models} {params.summary_stats}"

	
rule final:
	input:
		summary_fel="hyphy/{locus}/fel/{locus}.FEL.summary",
		summary_fubar="hyphy/{locus}/fubar/{locus}.FUBAR.summary",
		directory="codeml/{locus}/stats/"

	output:
		"hyphy/{locus}/done.txt"

	threads: 1

	shell:
		"touch {output}"
	
