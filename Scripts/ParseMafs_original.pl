#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;


$start = -1;
$stop = -1;
$scaff1 = "";
$scaff2 = "";
$length = 0;
$ave = 0;
 
while (<FILE>) {
	$line = $_; chomp($line);
	@tabs = split(/\t/,$line);


	if ($tabs[5] >= 65 && $start == -1) {
		$start = $tabs[1];
		$scaff1 = $tabs[0];
		$stop = $tabs[1];
		$scaff2 = $tabs[0];
		$length = 1;
		$ave = $tabs[5];
	} elsif ($tabs[5] >= 65 && $start != -1) {
		$stop = $tabs[1];
		$scaff2 = $tabs[0];
		$length++;
		$ave = $ave + $tabs[5];
	} elsif ($tabs[5] <= 65 && $start != -1) {
		$ave = $ave/$length;

		if ($length == 84 && $scaff1 eq $scaff2 && ($stop - $start) == 83) {
			print "$scaff1\t$start\t$scaff2\t$stop\t$length\t$ave\n";
		}

		$start = -1;
		$stop = -1;
		$scaff1 = "";
		$scaff2 = "";
		$length = 0;
		$ave = 0;
	}
}
close FILE;


