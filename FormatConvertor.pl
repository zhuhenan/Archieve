#!/usr/bin/perl -w
use strict;
use warnings;

# ------------------------------------------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------------------------------------------
# Import local tool library
use Getopt::Long;
use Pod::Usage;

# ------------------------------------------------------------------------------------------------------------------
# Main functions
# Call different method according to arguments
# ------------------------------------------------------------------------------------------------------------------

my $help = 0;
GetOptions(
	'help|h'     => \ $help,
	'method|m=s' => \ my $method,
	'gff3|g=s'   => \ my $gff,
	'rep|r=s',   => \ my $rep,
	'bed|b=s'    => \ my $bed
);
# Message for help function
pod2usage(-exitval => 0, -verbose => 99, -sections => [qw(NAME|SYNOPSIS|OPTIONS)]) if $help;

# Function controller
if ($method eq "Gff2Bed") {
    if (length($gff) && length($bed)) {
        &Gff2Bed($gff, $bed);
    } else {
		pod2usage(-exitval => 0, -verbose => 99, -sections => [qw(NAME|DESCRIPTION/Gff2Bed)]);
	}
} elsif ($method eq "RPout2Bed") {
	if (length($rep) && length($bed)) {
        &RPout2Bed($rep, $bed);
    } else {
		pod2usage(-exitval => 0, -verbose => 99, -sections => [qw(NAME|DESCRIPTION/RPout2Bed)]);
	}
} else {
	pod2usage(-exitval => 0, -verbose => 99, -sections => [qw(SYNOPSIS|OPTIONS)])
}

# ------------------------------------------------------------------------------------------------------------------
# Sub functions - Gff2Bed
# convert gff3 to bed
# ------------------------------------------------------------------------------------------------------------------
sub Gff2Bed {
	print "function in progress\n";
}

# ------------------------------------------------------------------------------------------------------------------
# Sub functions - RPout2Bed
# convert RepeatMasker output to bed
# ------------------------------------------------------------------------------------------------------------------
sub RPout2Bed {
	
	#
	my ($rep, $bed) = @_;
	
	#
	open(my $input,  "<", $rep) or die "$rep not exists";
	open(my $output, ">", $bed);
	#
	for (my $i = 0; $i < 3; $i ++) { $_ = <$input>; }
	
	#
	my $i = 1;
	while (defined(my $str = <$input>)) {
        $str =~ s/^\s+//;
		my @bed = split(" +", $str);
		#
		my $chrom  = $bed[4];
		my $start  = $bed[5];
		my $end    = $bed[6];
		my $name   = "Seq$i\_$bed[4]_$bed[5]_$bed[6]_$bed[8]_$bed[9]#$bed[10]";
		my $score  = $bed[0];
		my $strand = "";
		if ($bed[8] eq "+") {
            $strand = $bed[8];
        } else {
			$strand = "-";
		}
        
		#
		print $output "$chrom\t$start\t$end\t$name\t$score\t$strand\n";
		#
		$i ++;
    }
}


# ------------------------------------------------------------------------------------------------------------------
# Help message sections
# ------------------------------------------------------------------------------------------------------------------
__END__

=head1 NAME

	FormatConvertor - Henan's multiple formats convertor

=head1 SYNOPSIS

	FormatConvertor.pl -m <method> <options>

=head1 OPTIONS

	Gff2Bed          convert gff3 to bed
	RPout2Bed        convert RepeatMasker output to bed

=head1 DESCRIPTION

=head2 Gff2Bed

=item B<useage:>

	FormatConvertor.pl -m Gff2Bed -gff <gff3> -bed <bed>

=item B<options:>
	
	-gff|g          input gff3 file
	-bed|b          output bed file

=head2 RPout2Bed

=item B<useage:>

	FormatConvertor.pl -m RPout2Bed -rep <repeatmasker.out> -bed <bed>

=item B<options:>
	
	-rep|r          input RepeatMasker output file
	-bed|b          output bed file

=cut