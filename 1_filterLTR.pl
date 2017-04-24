use strict;

my @rps = glob("example/*.out");
foreach my $rp (@rps) {
	print $rp,"\n";
	
	open(my $input,  "<", $rp);
	open(my $output, ">", $rp.".LTR");
	
	for (my $i = 0; $i < 3; $i ++) {
		print $output $_;
	}
	
	while (defined($_ = <$input>)) {
        if ($_ =~ /LTR/ && $_ !~ /-int/) {
            print $output $_;
        }
    }
}