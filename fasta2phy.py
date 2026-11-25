from pathlib import Path

fasta = "Cervinae_10inds_1out.concat.fa"
output = "Cervinae_10inds_1out.concat.phy"

seqs = {}
name = None

with open(fasta) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seqs:
                raise ValueError(f"Duplicate sequence name: {name}")
            seqs[name] = ""
        else:
            seqs[name] += line

lengths = {len(s) for s in seqs.values()}

if len(lengths) != 1:
    raise ValueError("Sequences have different lengths!")

with open(output, "w") as f:
    f.write(f"{len(seqs)} {lengths.pop()}\n")
    for name, seq in seqs.items():
        f.write(f"{name:<20}{seq}\n")

print(f"Converted {len(seqs)} sequences to {output}")
