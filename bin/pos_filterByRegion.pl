# Author : SHujia Huang
# Modified: Dongbing Liu
# Date   : 2012-07-08, 2013-03-19
#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <target Region file(chr beg end ...)> <variant File(chr pos ...)>\n"if(@ARGV!=2);

my ( $TRFile, $VRFile) = @ARGV;

my %region;
my %index;

if ($TRFile =~ /\.gz$/) {open I, "gzip -dc $TRFile |" or die "Cannot open file : $TRFile\n";}
else {open I, $TRFile or die "Cannot open file : $TRFile\n";}
my (@tmp, $key);
while ( <I> ) {
	chomp;
	@tmp = split; # chr beg end ...
	$key = $tmp[0];
    my @tmp1=($tmp[1],$tmp[2]);
	push( @{ $region{ $key } }, [ @tmp1 ] );
}
close( I );

for my $key ( keys %region ) {

	@{$region{$key}} = sort { $a->[0] <=> $b->[0] } @{$region{$key}};
	$index{$key}     = 0;
}

##################
my $mapNum = 0;
if($VRFile=~/\.gz$/){open I,"gzip -dc $VRFile|" or die $!;}else{open I, $VRFile or die "Cannot open file : $VRFile\n";}
while ( <I> ) {
	chomp;
    my $line=$_;
    if($line=~/^#/){print "$line\n";next;}
	@tmp = split /\s+/,$line; # chr pos ... , sorted by coordinate
	next if ( !defined $region{$tmp[0]} );
    
    my $flag=1;
    my $tmp=$tmp[2];
    $tmp[2]=$tmp[1];
	for ( my $i = $index{$tmp[0]}; $i < @{ $region{$tmp[0]} }; ++$i  ) {

		last if $tmp[2] < $region{$tmp[0]}[$i][0];
		next if $tmp[1] > $region{$tmp[0]}[$i][1];

		if ($flag){
			$index{$tmp[0]} = $i;
			$flag = 0;
		}

		print "$line\n";
	}
}
close ( I );


