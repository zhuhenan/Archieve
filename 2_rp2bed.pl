use strict;

my @rps = glob("example/*.out.LTR");
foreach my $rp (@rps) {
	print $rp,"\n";
	$_ = $rp;
	s/repeatmasker/beds/;
	system "perl FormatConvertor.pl -m RPout2Bed -r $rp -b $_.bed";
}