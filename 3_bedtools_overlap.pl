use strict;

my @beds = glob("example/*.LTR.bed");
foreach my $bed (@beds) {
	print $bed,"\n";
	system "sort -k1,1 -k2,2n $bed > $bed.sort";
	
	system "bedtools merge -i $bed.sort -c 1,4 -o count,collapse > $bed.sort.merge";
}