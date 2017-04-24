use strict;

my $genome_dir = "~/Documents/application/DIGS-tool/genome/Mammalia/Equus_caballus/complete/goldenpath_equCab2_Dec13";
my ($input, $output, %seq, $group);

my @seqs = glob("example/*.sort.solo.fa");
foreach my $file (@seqs) {
	print $file,"\n";
	open($input, "<", $file);
	while (defined($_ = <$input>)) {
        if ($_ =~ />/) {
            @_ = split("_", $_);
			$group = $2 if ($_ =~ /(C|\+)_(.+)#LTR/);
			$seq{$group} .= ">$_[1]_$_[2]_$_[3]_$group\n"
        } else {
			$seq{$group} .= "$_\n";
		}
    }
}

foreach my $key (keys %seq) {
	print $key,"\n";
	open($output, ">", "example/$key.fa");
	print $output $seq{$key};
}