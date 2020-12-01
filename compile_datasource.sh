#!/usr/bin/env bash
# Jbrowse2 datasource deployment script
# Author Muhammad Elhossary | elhossary@zbmed.de
#set -e
WORKING_DIR="$@"
main(){
	## Setting variables
	BIN_DIR=bin
	RAW_DATA_DIR=${WORKING_DIR}
	REFSEQ_DIR=${RAW_DATA_DIR}/reference_sequence
	ANNOTATIONS_DIR=${RAW_DATA_DIR}/annotations
	COVERAGE_DIR=${RAW_DATA_DIR}/coverage
	DATA_POOL=data_pool
	CHROM_SIZES_FILE="${DATA_POOL}/all.chrom.sizes"
	BW_LIST=*.bw
	ASSEMBLY_NAME=$(basename ${WORKING_DIR})
	ASSEMBLY_ALIAS=${ASSEMBLY_NAME}
	mkdir ${BIN_DIR}
	mkdir ${DATA_POOL}

	## Function calls
	index_fasta_file
	add_assembly_and_refseq_track
	split_gff_to_forward_reverse
	add_annotation_tracks
	#download_wigToBigWig_tool
	#convert_wig_to_bigwig
	add_coverage_track
}
index_fasta_file(){
	for file_name in "${REFSEQ_DIR}"/*.fa
	do
		echo "==> Indexing fasta file: $(basename ${file_name} .fa)"
		samtools faidx ${file_name}
	done
}
add_assembly_and_refseq_track(){
	for file_name in "${REFSEQ_DIR}"/*.fa
	do
		echo "==> Creating assembly: $(basename ${file_name} .fa)"
		jbrowse add-assembly ${file_name} \
		--load copy \
		--type indexedFasta \
		--name ${ASSEMBLY_NAME} \
		--alias ${ASSEMBLY_ALIAS} \
		--force \
		--out ${DATA_POOL}
	done
}
split_gff_to_forward_reverse(){
	echo "==> Splitting GFF files for: ${ORGANISM_DATASOURCE_DIR}"
	for file_name in "${ANNOTATIONS_DIR}"/*.gff
	do
		if ([[ ! $file_name = *"reverse"* ]] && [[ ! $file_name = *"forward"* ]])
		then
			echo -e "===> Splitting annotation: $(basename ${file_name})"
			grep -P "\t-\t" ${file_name} > "$(dirname ${file_name})/$(basename ${file_name} .gff)_reverse.gff"
			grep -Pv "\t-\t" ${file_name} > "$(dirname ${file_name})/$(basename ${file_name} .gff)_forward.gff"
			rm "${file_name}"
		fi
	done
}
add_annotation_tracks(){
	for FILE in "${ANNOTATIONS_DIR}"/*.gff
	do
		bname=$(basename ${FILE} .gff)
		dname=$(dirname ${FILE})
		# Index
		echo "Indexing track ${bname}"
		(grep ^"#" ${FILE}; grep -v ^"#" ${FILE} | sort -k1,1 -k4,4n) | bgzip > "${dname}/${bname}.gff.gz";
		tabix -p gff "${dname}/${bname}.gff.gz";
		# Add
		jbrowse add-track "${dname}/${bname}.gff.gz" \
		--assemblyNames ${ASSEMBLY_NAME} \
		--description "Annotations" \
		--force \
		--load copy \
		--name ${bname} \
		--category "Annotations" \
		--out ${DATA_POOL}
	done
}
download_wigToBigWig_tool(){
	echo "Downloading tool wigToBigWig"
	wget -P "${BIN_DIR}" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
	chmod u+x "${BIN_DIR}"/wigToBigWig

}
update_chromosomes_sizes_file(){
	echo "==> Updating chromosomes sizes file for: ${ORGANISM_DATASOURCE_DIR}"
	for file_name in "${REFSEQ_DIR}"/*.*
	do
		# Get all chromosomes sizes
		cat "${file_name}" | cut -d' ' -f 1 | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >> "${CHROM_SIZES_FILE}"
		# Remove duplicates if found
		sort < "${CHROM_SIZES_FILE}" | uniq > "${CHROM_SIZES_FILE}~" && mv "${CHROM_SIZES_FILE}~" "${CHROM_SIZES_FILE}"
	done
}
convert_wig_to_bigwig(){
	echo "==> Converting wig to bigwig for: ${ORGANISM_DATASOURCE_DIR}"
	update_chromosomes_sizes_file
	for FILE in "${COVERAGE_DIR}"/*.wig
	do
		echo -e "===> Converting: $FILE"
		"${BIN_DIR}"/wigToBigWig "$FILE" "${CHROM_SIZES_FILE}" "${COVERAGE_DIR}"/$(basename "$FILE" .wig).bw
		rm "${FILE}"
	done
}
add_coverage_track(){
	for FILE in "${COVERAGE_DIR}"/*.bw
	do
		jbrowse add-track ${FILE} \
		--assemblyNames ${ASSEMBLY_NAME} \
		--description "Coverage" \
		--force \
		--load copy \
		--name $(basename "$FILE" .bw) \
		--category "Coverage" \
		--out ${DATA_POOL}
	done
}

main
