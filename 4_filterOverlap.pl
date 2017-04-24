use strict;

my ($input, $output, %hash);

my @sorts = glob("example/*.sort");
foreach my $sort (@sorts) {
	#
	print $sort,"\n";
	open($input, "<", $sort);
	#
	undef %hash;
	while (defined($_ = <$input>)) {
        @_ = split("\t", $_);
		$hash{$_[3]} = $_[4];
    }
    #
	open($input,  "<", "$sort.merge");
	open($output, ">", "$sort.solo");
	#
	while (defined($_ = <$input>)) {
		chomp $_;
        @_ = split("\t", $_);
		#
		print $output "$_[0]\t$_[1]\t$_[2]\t";
		#
		if ($_[3] > 1) {
            @_ = split(",", $_[-1]);
			#
			my $best = 0;
			my $label = "";
			foreach my $alt (@_) {
				if ($best < $hash{$alt}) {
                    $best  = $hash{$alt};
					$label = $alt;
                }
			}
			#
			print $output $label."\t".&strand($label)."\n";
			#
        } else {
			print $output "$_[-1]\t".&strand($_[-1])."\n";
		}
    }
}

sub strand {
	my $label = shift;
	my $strand = "+";
	if ($label !~ /\+/) {
        $strand = "-";
    }
	return $strand;
}