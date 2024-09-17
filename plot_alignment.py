#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

# Alignment stats file and output file
stats_file = sys.argv[1]
output_file = sys.argv[2]

# Parse alignment statistics
with open(stats_file) as f:
    lines = f.readlines()

total_reads = int(lines[0].split()[0])
mapped_reads = int(lines[4].split()[0])
unmapped_reads = total_reads - mapped_reads

# Create a DataFrame for visualization
data = {
    'Metric': ['Total reads', 'Mapped reads', 'Unmapped reads'],
    'Value': [total_reads, mapped_reads, unmapped_reads]
}
df = pd.DataFrame(data)
df['Percentage'] = df['Value'] / df['Value'].sum() * 100

# Create a barplot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
plot = sns.barplot(x="Metric", y="Percentage", data=df, palette="Blues_d")
plt.title("Sequenz-Übereinstimmung mit der Referenz")
plt.ylabel("Prozent (%)")
plt.xlabel("Metrik")

# Add percentages on top of bars
for index, row in df.iterrows():
    plot.text(index, row['Percentage'] + 1, f'{row["Percentage"]:.2f}%', color='black', ha="center")

# Save the plot as a PNG image
plt.savefig(output_file)
plt.show()

