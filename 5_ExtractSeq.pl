use strict;

my $genome_dir = "~/Documents/application/DIGS-tool/genome/Mammalia/Equus_caballus/complete/goldenpath_equCab2_Dec13";
my ($input, $output, %seq, $group);

my @solos = glob("example/*.sort.solo");
foreach my $solo (@solos) {
	print $solo,"\n";
	my $chrom = $1 if ($solo =~ /(chr\w+)\.fa/);
	my $in_fa = $genome_dir."/$chrom.fa";
	print $in_fa,"\n";
	
	my $command = "bedtools getfasta -fi $in_fa -bed $solo -fo $solo.fa -name -s";
	system "$command";
}