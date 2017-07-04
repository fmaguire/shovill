my $MIN_K = 31;   # sensible minimum, although prefer higher
my $MAX_K = 127;  # spades max
my $MINBQ = 3;    # for trimming and stats
my $R1 = $ARGV[0];

sub read_file {
    my $fn = shift;
    open my $fh, '<', $fn or die "Can't open file $!";
    chomp(my @lines = <$fh>);
    return @lines;
}

sub read_stats {
  my($R1) = @_;
  my $outfile = "10-seqtk.tab";
  my $stat;
  run_cmd("seqtk fqchk -q $MINBQ \Q$R1\E > $outfile");
  my @row = read_file($outfile);
#  msg($row[0]);
  for my $tag ('min_len', 'max_len', 'avg_len') {
    $row[0] =~ m/$tag:\s*(\d+(\.\d+)?);/ or err("Can't parse $tag from: $row[0]");
    $stat->{$tag} = int( $1 + 0.5 );
  }
  $row[2] =~ m/^ALL\s+(\d+)/ or err("Can't parse ALL #bases from: $row[2]");
  $stat->{'total_bp'} = $1 * 2;  # multiply by 2 as only using R1 (hack)
  return $stat;
}

sub err {
  msg(@_);
  exit(1);
}

sub msg {
  my $msg = "@_\n";
  print STDERR $msg;
  push @LOG, $msg;
}

sub run_cmd {
  my($cmd, $logfile) = @_;
  if ($logfile) {
    $cmd .= $debug ? " 2>&1 | tee -a $logfile" : " >> $logfile 2>&1";
  }
  msg("Running: $cmd");
  system($cmd)==0 or err("Error $? running command");
}

sub main {
    my $cpus = 2;
    
    msg("Calculating read stats");
    my $stat = read_stats($R1);

    msg("Calculating gsize");
    my $minkc = 3;  
    run_cmd("kmc -ci$minkc -k25 -t$cpus \Q$R1\E kmc /tmp", "20-kmc.log");
    my($ucount) = grep { m/unique\s+counted/ } read_file('20-kmc.log');
    $ucount =~ m/(\d+)$/ or err("Could not determine unique counted k-mers using kmc");

    $gsize = $1;

    open(my $fh, '>', "genome_size.txt");
    print $fh $gsize;
    close $fh

   
    my $depth = int( $stat->{'total_bp'} / $gsize );
    
    my $RLEN = $stat->{'avg_len'};
    msg("Average read length looks like $RLEN bp");
    $MAX_K = min( $MAX_K, int(0.8 * $RLEN) );  
    msg("Setting maximum k-mer to $MAX_K");
    my $minlen = $RLEN;
    
    # Choosing some kmers
    $MIN_K = 21 if $stat->{'avg_len'} < 75;
    my $kn = 5; 
    my $ks = max(5, int( ($MAX_K - $MIN_K) / ($kn-1) ) );
    $ks++ if $ks % 2 == 1; # need even step to ensure odd values only
    my @kmers;
    for (my $k=$MIN_K; $k <= $MAX_K; $k+=$ks) {
        push @kmers, $k;
    }
    msg("Estimated K-mers: @kmers [kn=$kn, ks=$ks, kmin=$MIN_K, kmax=$MAX_K]");
    $kmers = join(',', @kmers);

    open(my $fh, '>', "kmers.txt");
    print $fh $kmers;
    close $fh

    msg("Using kmers: $kmers");
}

&main();
