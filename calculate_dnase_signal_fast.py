promoter_file = "sorted_chr_promoters_with_ID.bed"
dnase_file = "wgEncodeUwDnaseGm12878RawRep1.bigWig.bgr"
output_file = "promoter_dnase_signal.txt"

dnase = []

with open(dnase_file) as f:
    for line in f:
        chrom, start, end, signal = line.strip().split("\t")
        dnase.append((chrom, int(start), int(end), float(signal)))

dnase_index = 0

with open(promoter_file) as p, open(output_file, "w") as out:

    for line in p:
        gene, chrom, start, end = line.strip().split("\t")
        start = int(start)
        end = int(end)

        total_signal = 0

        while dnase_index < len(dnase):
            dchrom, dstart, dend, signal = dnase[dnase_index]

           
            if dend <= start:
                dnase_index += 1
                continue

            
            if dstart >= end:
                break

            
            overlap_start = max(start, dstart)
            overlap_end = min(end, dend)

            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                total_signal += overlap_length * signal

            dnase_index += 1

        out.write(f"{gene}\t{total_signal}\n")

print("Finished!")
