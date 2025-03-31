def clean_sequence(sequence)
  sequence.gsub(/[^atcgATCG]/, '').downcase 
end

# Функция для поиска минимальной специфичной последовательности
def find_min_unique(seq1, seq2)
  (1..seq1.length).each do |length|
    (0..seq1.length - length).each do |i|
      substring = seq1[i, length]
      return substring unless seq2.include?(substring)
    end
  end
  nil
end

# Функция для поиска максимальной общей подпоследовательности (LCS)
def find_lcs(seq1, seq2)
  m = seq1.length
  n = seq2.length
  lcs_table = Array.new(m + 1) { Array.new(n + 1, 0) }

  (1..m).each do |i|
    (1..n).each do |j|
      if seq1[i - 1] == seq2[j - 1]
        lcs_table[i][j] = lcs_table[i - 1][j - 1] + 1
      else
        lcs_table[i][j] = [lcs_table[i - 1][j], lcs_table[i][j - 1]].max
      end
    end
  end

  # Восстанавливаем LCS из таблицы
  lcs_length = lcs_table[m][n]
  lcs = ''
  i, j = m, n

  while i > 0 && j > 0
    if seq1[i - 1] == seq2[j - 1]
      lcs = seq1[i - 1] + lcs
      i -= 1
      j -= 1
    elsif lcs_table[i - 1][j] > lcs_table[i][j - 1]
      i -= 1
    else
      j -= 1
    end
  end

  [lcs, lcs_length.to_f / [m, n].max] # Возвращаем LCS и отношение длины к общей длине генома
end

covid_sequence = clean_sequence(<<~HEREDOC).slice(0, 6000)
attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
 gttctctaaa cgaactttaa aatctgtgtg gctgtcactc ggctgcatgc ttagtgcact
 cacgcagtat aattaataac taattactgt cgttgacagg acacgagtaa ctcgtctatc
 ttctgcaggc tgcttacggt ttcgtccgtg ttgcagccga tcatcagcac atctaggttt
 cgtccgggtg tgaccgaaag gtaagatgga gagccttgtc cctggtttca acgagaaaac
 acacgtccaa ctcagtttgc ctgttttaca ggttcgcgac gtgctcgtac gtggctttgg
 agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
 cttagtagaa gttgaaaaag gcgttttgcc tcaacttgaa cagccctatg tgttcatcaa
 acgttcggat gctcgaactg cacctcatgg tcatgttatg gttgagctgg tagcagaact
 cgaaggcatt cagtacggtc gtagtggtga gacacttggt gtccttgtcc ctcatgtggg
 cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
 tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
 tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg ttacccgtga
 actcatgcgt gagcttaacg gaggggcata cactcgctat gtcgataaca acttctgtgg
 ccctgatggc taccctcttg agtgcattaa agaccttcta gcacgtgctg gtaaagcttc
 atgcactttg tccgaacaac tggactttat tgacactaag aggggtgtat actgctgccg
 tgaacatgag catgaaattg cttggtacac ggaacgttct gaaaagagct atgaattgca
 gacacctttt gaaattaaat tggcaaagaa atttgacacc ttcaatgggg aatgtccaaa
 ttttgtattt cccttaaatt ccataatcaa gactattcaa ccaagggttg aaaagaaaaa
 gcttgatggc tttatgggta gaattcgatc tgtctatcca gttgcgtcac caaatgaatg
 caaccaaatg tgcctttcaa ctctcatgaa gtgtgatcat tgtggtgaaa cttcatggca
 gacgggcgat tttgttaaag ccacttgcga attttgtggc actgagaatt tgactaaaga
 aggtgccact acttgtggtt acttacccca aaatgctgtt gttaaaattt attgtccagc
 atgtcacaat tcagaagtag gacctgagca tagtcttgcc gaataccata atgaatctgg
 cttgaaaacc attcttcgta agggtggtcg cactattgcc tttggaggct gtgtgttctc
 ttatgttggt tgccataaca agtgtgccta ttgggttcca cgtgctagcg ctaacatagg
 ttgtaaccat acaggtgttg ttggagaagg ttccgaaggt cttaatgaca accttcttga
 aatactccaa aaagagaaag tcaacatcaa tattgttggt gactttaaac ttaatgaaga
 gatcgccatt attttggcat ctttttctgc ttccacaagt gcttttgtgg aaactgtgaa
 aggtttggat tataaagcat tcaaacaaat tgttgaatcc tgtggtaatt ttaaagttac
 aaaaggaaaa gctaaaaaag gtgcctggaa tattggtgaa cagaaatcaa tactgagtcc
 tctttatgca tttgcatcag aggctgctcg tgttgtacga tcaattttct cccgcactct
 tgaaactgct caaaattctg tgcgtgtttt acagaaggcc gctataacaa tactagatgg
 aatttcacag tattcactga gactcattga tgctatgatg ttcacatctg atttggctac
 taacaatcta gttgtaatgg cctacattac aggtggtgtt gttcagttga cttcgcagtg
 gctaactaac atctttggca ctgtttatga aaaactcaaa cccgtccttg attggcttga
 agagaagttt aaggaaggtg tagagtttct tagagacggt tgggaaattg ttaaatttat
 ctcaacctgt gcttgtgaaa ttgtcggtgg acaaattgtc acctgtgcaa aggaaattaa
 ggagagtgtt cagacattct ttaagcttgt aaataaattt ttggctttgt gtgctgactc
 tatcattatt ggtggagcta aacttaaagc cttgaattta ggtgaaacat ttgtcacgca
 ctcaaaggga ttgtacagaa agtgtgttaa gttaaatcca gagaagaaac tggcctactc
 atgcc
 tctaaaagcc ccaaaagaaa ttatcttctt agagggagaa acacttccca cagaagtgtt
 aacagaggaa gttgtcttga aaactggtga tttacaacca ttagaacaac ctactagtga
 agctgttgaa gctccattgg ttggtacacc agtttgtatt aacgggctta tgttgctcga
 aatcaaagac acagaaaagt actgtgccct tgcacctaat atgatggtaa caaacaatac
 cttcacactc aaaggcggtg caccaacaaa ggttactttt ggtgatgaca ctgtgataga
 agtgcaaggt tacaagagtg tgaatatcac ttttgaactt gatgaaagga ttgataaagt
 acttaatgag aagtgctctg cctatacagt tgaactcggt acagaagtaa atgagttcgc
 ctgtgttgtg gcagatgctg tcataaaaac tttgcaacca gtatctgaat tacttacacc
 actgggcatt gatttagatg agtggagtat ggctacatac tacttatttg atgagtctgg
 tgagtttaaa ttggcttcac atatgtattg ttctttctac cctccagatg aggatgaaga
 agaaggtgat tgtgaagaag aagagtttga gccatcaact caatatgagt atggtactga
 agatgattac caaggtaaac ctttggaatt tggtgccact tctgctgctc ttcaacctga
 agaagagcaa gaagaagatt ggttagatga tgatagtcaa caaactgttg gtcaacaaga
 cggcagtgag gacaatcaga caactactat tcaaacaatt gttgaggttc aacctcaatt
 agagatggaa cttacaccag ttgttcagac tattgaagtg aatagtttta gtggttattt
 aaaacttact gacaatgtat acattaaaaa tgcagacatt gtggaagaag ctaaaaaggt
 aaaaccaaca gtggttgtta atgcagccaa tgtttacctt aaacatggag gaggtgttgc
 aggagcctta aataaggcta ctaacaatgc catgcaagtt gaatctgatg attacatagc
 tactaatgga ccacttaaag tgggtggtag ttgtgtttta agcggacaca atcttgctaa
 acactgtctt catgttgtcg gcccaaatgt taacaaaggt gaagacattc aacttcttaa
 gagtgcttat gaaaatttta atcagcacga agttctactt gcaccattat tatcagctgg
 tatttttggt gctgacccta tacattcttt aagagtttgt gtagatactg ttcgcacaaa
 tgtctactta gctgtctttg ataaaaatct ctatgacaaa cttgtttcaa gctttttgga
 aatgaagagt gaaaagcaag ttgaacaaaa gatcgctgag attcctaaag aggaagttaa
 gccatttata actgaaagta aaccttcagt tgaacagaga aaacaagatg ataagaaaat
 caaagcttgt gttgaagaag ttacaacaac tctggaagaa actaagttcc tcacagaaaa
 cttgttactt tatattgaca ttaatggcaa tcttcatcca gattctgcca ctcttgttag
 tgacattgac atcactttct taaagaaaga tgctccatat atagtgggtg atgttgttca
 agagggtgtt ttaactgctg tggttatacc tactaaaaag gctggtggca ctactgaaat
 gctagcgaaa gctttgagaa aagtgccaac agacaattat ataaccactt acccgggtca
 gggtttaaat ggttacactg tagaggaggc aaagacagtg cttaaaaagt gtaaaagtgc
 cttttacatt ctaccatcta ttatctctaa tgagaagcaa gaaattcttg gaactgtttc
 ttggaatttg cgagaaatgc ttgcacatgc agaagaaaca cgcaaattaa tgcctgtctg
 tgtggaaact aaagccatag tttcaactat acagcgtaaa tataagggta ttaaaataca
 agagggtgtg gttgattatg gtgctagatt ttacttttac accagtaaaa caactgtagc
 gtcacttatc aacacactta acgatctaaa tgaaactctt gttacaatgc cacttggcta
 tgtaacacat ggcttaaatt tggaagaagc tgctcggtat atgagatctc tcaaagtgcc
 agctacagtt tctgtttctt cacctgatgc tgttacagcg tataatggtt atcttacttc
 ttcttctaaa acacctgaag aacattttat tgaaaccatc tcacttgctg gttcctataa
 agattggtcc tattctggac aatctacaca actaggtata gaatttctta agagaggtga
 taaaagtgta tattacacta gtaatcctac cacattccac ctagatggtg aagttatcac
 ctttgacaat cttaagacac ttctttcttt gagagaagtg aggactatta aggtgtttac
 aacagtagac aacattaacc tccacacgca agttgtggac atgtcaatga catatggaca
 acagtttggt ccaacttatt tggatggagc tgatgttact aaaataaaac ctcataattc
HEREDOC

h5n1_sequence = clean_sequence(<<~HEREDOC).slice(0, 6000)
agcaaaagca gggtagataa tcactcaatg agtgacatcg aagccatggc atctcaaggc
 accaaacgat catatgaaca aatggagact ggtggggaac gccaggatgc cacagaaatc
 agagcatctg tcggaaaaat gattggtgga atcgggaaat tctacatcca aatgtgcact
 gaactcaaac tcagtgacta tgagggacga ctaatccaaa atagcataac aatagagaga
 atggtgctct ctgcttttga tgagagaaga aataaatacc tagaagagca tcccagtgct
 gggaaggatc ctaagaaaac tggaggaccc atatatagaa gagtagacgg aaagtggatg
 agagaactca ttctttatga caaagaagaa ataaggagag tttggcgcca agcaaataat
 ggtgaagatg caacagctgg tcttactcat atcatgattt ggcattccaa tctgaatgat
 gccacgtacc agagaacaag agcgcttgtt cgcaccggaa tggatcccag aatgtgctct
 ctaatgcaag gttcaacact tcccagaagg tctggggccg caggtgctgc agtgaaagga
 gttggaacaa tagcaatgga attgatcaga atgatcaaac gtgggatcaa tgaccgaaac
 ttctggagag gtgaaaatgg acgaaggaca agggttgcat atgaaagaat gtgcaatatt
 ctcaaaggaa aatttcagac agctgcccag agggcaatga tggatcaagt gagagaaagt
 cggaacccag gaaacgctga gattgaagat ctcattttcc tggcacggtc agcacttatt
 ctaaggggat cagttgcaca taagtcttgc ctgcctgctt gtgtgtatgg gcttgcagtg
 gcaagtgggc atgactttga aagggaaggg tattcactgg tcgggataga cccgtttaaa
 ttactccaaa acagtcaagt gttcagcctg ataagaccaa atgaaaaccc agctcacaag
 agtcaattag tgtggatggc atgccactct gctgcatttg aggatctaag ggtatcaagt
 ttcataagag ggaagaaagt gattccaaga ggaaagcttt ccacaagagg agttcagatt
 gcttcaaatg agaatgtgga agccatggat tccaatacct tagagctgag aagcagatac
 tgggccataa ggaccagaag tggaggaaat accaatcaac aaaaagcatc cgctggccag
 atcagtgtgc aacctacatt ctcagtgcaa cggaatctcc cttttgaaag agcaaccgtt
 atggcagcct tcagcgggaa caatgaagga cggacatccg atatgcgaac agaagttata
 aagatgatgg aaagtgcaaa gccagaagat ttgtccttcc aggggcgggg agtcttcgag
 ctctcggacg aaaaggcaac gagcccgatc gtgccttcct ttgacatgag taatgaagga
 tcttatttct tcggagacaa tgcagaggag tatgacagtt gaggaaaaat acccttgttt
 ctact
HEREDOC
delta_sequence = clean_sequence(<<~HEREDOC).slice(0, 6000)
caaaccaacc aactttcgat ctcttgtaga tctgttctct aaacgaactt taaaatctgt
 gtggctgtca ctcggctgca tgcttagtgc actcacgcag tataattaat aactaattac
 tgtcgttgac aggacacgag taactcgtct atcttctgca ggctgcttac ggtttcgtcc
 gttttgcagc cgatcatcag cacatctagg ttttgtccgg gtgtgaccga aaggtaagat
 ggagagcctt gtccctggtt tcaacgagaa aacacacgtc caactcagtt tgcctgtttt
 acaggttcgc gacgtgctcg tacgtggctt tggagactcc gtggaggagg tcttatcaga
 ggcacgtcaa catcttaaag atggcacttg tggcttagta gaagttgaaa aaggcgtttt
 gcctcaactt gaacagccct atgtgttcat caaacgttcg gatgctcgaa ctgcacctca
 tggtcatgtt atggttgagc tggtagcaga actcgaaggc attcagtacg gtcgtagtgg
 tgagacactt ggtgtccttg tccctcatgt gggcgaaata ccagtggctt accgcaaggt
 tcttcttcgt aagaacggta ataaaggagc tggtggccat agttacggcg ccgatctaaa
 gtcatttgac ttaggcgacg agcttggcac tgatccttat gaagattttc aagaaaactg
 gaacactaaa catagcagtg gtgttacccg tgaactcatg cgtgagctta acggaggggc
 atacactcgc tatgtcgata acaacttctg tggccctgat ggctaccctc ttgagtgcat
 taaagacctt ctagcacgtg ctggtaaagc ttcatgcact ttgtccgaac aactggactt
 tattgacact aagaggggtg tatactgctg ccgtgaacat gagcatgaaa ttgcttggta
 cacggaacgt tctgaaaaga gctatgaatt gcagacacct tttgaaatta aattggcaaa
 gaaatttgac accttcaatg gggaatgtcc aaattttgta tttcccttaa attccataat
 caagactatt caaccaaggg ttgaaaagaa aaagcttgat ggctttatgg gtagaattcg
 atctgtctat ccagttgcgt caccaaatga atgcaaccaa atgtgccttt caactctcat
 gaagtgtgat cattgtggtg aaacttcatg gcagacgggc gattttgtta aagccacttg
 cgaattttgt ggcactgaga atttgactaa agaaggtgcc actacttgtg gttacttacc
 ccaaaatgct gttgttaaaa tttattgtcc agcatgtcac aattcagaag taggacctga
 gcatagtctt gccgaatacc ataatgaatc tggcttgaaa accattcttc gtaagggtgg
 tcgcactatt gcctttggag gctgtgtgtt ctcttatgtt ggttgccata acaagtgtgc
 ctattgggtt ccacgtgcta gcgctaacat aggttgtaac catacaggtg ttgttggaga
 aggttccgaa ggtcttaatg acaaccttct tgaaatactc caaaaagaga aagtcaacat
 caatattgtt ggtgacttta aacttaatga agagatcgcc attattttgg catctttttc
 tgcttccaca agtgcttttg tggaaactgt gaaaggtttg gattataaag cattcaaaca
 aattgttgaa tcctgtggta attttaaagt tacaaaagga aaagctaaaa aaggtgcctg
 gaatattggt gaacagaaat caatactgag tcctctttat gcatttgcat cagaggctgc
 tcgtgttgta cgatcaattt tctcccgcac tcttgaaact gctcaaaatt ctgtgcgtgt
 tttacagaag gccgctataa caatactaga tggaatttca cagtattcac tgagactcat
 tgatgctatg atgttcacat ctgatttggc tactaacaat ctagttgtaa tggcctacat
 tacaggtggt gttgttcagt tgacttcgca gtggctaact aacatctttg gcactgttta
 tgaaaaactc aaacccgtcc ttgattggct tgaagagaag tttaaggaag gtgtagagtt
 tcttagagac ggttgggaaa ttgttaaatt tatctcaacc tgtgcttgtg aaattgtcgg
 tggacaaatt gtcacctgtg caaaggaaat taaggagagt gttcagacat tctttaagct
 tgtaaataaa tttttggctt tgtgtgctga ctctatcatt attggtggag ctaaacttaa
 agccttgaat ttaggtgaaa catttgtcac gcactcaaag ggattgtaca gaaagtgtgt
 taaatccaga gaagaaactg gcctactcat gcctctaaaa gccccaaaag aaattatctt
 cttagaggga gaaacacttc ccacagaagt gttaacagag gaagttgtct tgaaaactgg
 tgatttacaa ccattagaac aacctactag tgaagctgtt gaagctccat tggttggtac
 accagtttgt attaacgggc ttatgttgct cgaaatcaaa gacacagaaa agtactgtgc
 ccttgcacct aatatgatgg taacaaacaa taccttcaca ctcaaaggcg gtgcaccaac
 aaaggttact tttggtgatg acactgtgat agaagtgcaa ggttacaaga gtgtgaatat
 cacttttgaa cttgatgaaa ggattgataa agtacttaat gagaagtgct ctgcctatac
 agttgaactc ggtacagaag taaatgagtt cgcctgtgtt gtggcagatg ctgtcataaa
 aactttgcaa ccagtatctg aattacttac accactgggc attgatttag atgagtggag
 tatggctaca tactacttat ttgatgagtc tggtgagttt aaattggctt cacatatgta
 ttgttctttt taccctccag atgaggatga agaagaaggt gattgtgaag aagaagagtt
 tgagccatca actcaatatg agtatggtac tgaagatgat taccaaggta aacctttgga
 atttggtgcc acttctgctg ctcttcaacc tgaagaagag caagaagaag attggttaga
 tgatgatagt caacaaactg ttggtcaaca agacggcagt gaggacaatc agacaactac
 tattcaaaca attgttgagg ttcaacctca attagagatg gaacttacac cagttgttca
 gactattgaa gtgaatagtt ttagtggtta tttaaaactt actgacaatg tatacattaa
 aaatgcagac attgtggaag aagctaaaaa ggtaaaacca acagtggttg ttaatgcagc
 caatgtttac cttaaacatg gaggaggtgt tgcaggagcc ttaaataagg ctactaacaa
 tgccatgcaa gttgaatctg atgattacat agctactaat ggaccactta aagtgggtgg
 tagttgtgtt ttaagcggac acaatcttgc taaacactgt cttcatgttg tcggcccaaa
 tgttaacaaa ggtgaagaca ttcaacttct taagagtgct tatgaaaatt ttaatcagca
 cgaagttcta cttgcaccat tattatcagc tggtattttt ggtgctgacc ctatacattc
 tttaagagtt tgtgtagata ctgttcgcac aaatgtctac ttagctgtct ttgataaaaa
 tctctatgac aaacttgttt caagcttttt ggaaatgaag agtgaaaagc aagttgaaca
 aaagatcgct gagattccta aagaggaagt taagccattt ataactgaaa gtaaaccttc
 agttgaacag agaaaacaag atgataagaa aatcaaagct tgtgttgaag aagttacaac
 aactctggaa gaaactaagt tcctcacaga aaacttgtta ctttatattg acattaatgg
 caatcttcat ccagattctg ccactcttgt tagtgacatt gacatcactt tcttaaagaa
 agatgctcca tatatagtgg gtgatgttgt tcaagagggt gttttaactg ctgtggttat
 acctactaaa aagtctggtg gcactactga aatgctagcg aaagctttga gaaaagtgcc
 aacagacaat tatataacca cttacccggg tcagggttta aatggttaca ctgtagagga
 ggcaaagaca gtgcttaaaa agtgtaaaag tgccttttac attctaccat ctattatctc
 taatgagaag caagaaattc ttggaactgt ttcttggaat ttgcgagaaa tgcttgcaca
 tgcagaagaa acacgcaaat taatgcctgt ctgtgtggaa actaaagcca tagtttcaac
 tatacagcgt aaatataagg gtattaaaat acaagagggt gtggttgatt atggtgctag
 attttacttt tacaccagta aaacaactgt agcgtcactt atcaacacac ttaacgatct
 aaatgaaact cttgttacaa tgccacttgg ctatgtaaca catggcttaa atttggaaga
 agctgctcgg tatatgagat ctctcaaagt gccagctaca gtttctgttt cttcacctga
 tgctgttaca gcgtataatg gttatcttac ttcttcttct aaaacacctg aagaacattt
 tattgaaacc atctcacttg ctggttccta taaagattgg tcctattctg gacaatctac
 acaactaggt atagaatttc ttaagagagg tgataaaagt gtatattaca ctagtaatcc
 taccacattc cacctagatg gtgaagttat cacctttgac aatcttaaga cacttctttc
 tttgagagaa gtgaggacta ttaaggtgtt tacaacagta gacaacatta acctccacac
 gcaagttgtg gacatgtcaa tgacatatgg acaacagttt ggtccaactt atttggatgg
 agctgatgtt actaaaataa aacctcataa ttcacatgaa ggtaaaacat tttatgtttt
HEREDOC

covid_cleaned = clean_sequence(covid_sequence)
h5n1_cleaned = clean_sequence(h5n1_sequence)
delta_cleaned = clean_sequence(delta_sequence)

min_unique_covid_h5n1 = find_min_unique(covid_cleaned, h5n1_cleaned)
min_unique_h5n1_covid = find_min_unique(h5n1_cleaned, covid_cleaned)

min_unique_covid_h5n1_delta = find_min_unique(covid_cleaned, h5n1_cleaned + delta_cleaned)
min_unique_delta_covid_h5n1 = find_min_unique(delta_cleaned, h5n1_cleaned + covid_cleaned)

max_common_covid_delta, ratio_lcs_covid_delta = find_lcs(covid_cleaned, delta_cleaned)

puts "Минимальная специфичная последовательность для COVID-19 (Ухань), отсутствующая в H5N1: #{min_unique_covid_h5n1}"
puts "Минимальная специфичная последовательность для H5N1, отсутствующая в COVID-19 (Ухань): #{min_unique_h5n1_covid}"
puts "Минимальная специфичная последовательность для COVID-19 (Ухань), отсутствующая в Delta и H5N1: #{min_unique_covid_h5n1_delta}"
puts "Минимальная специфичная последовательность для Delta, отсутствующая в COVID-19 (Ухань) и H5N1: #{min_unique_delta_covid_h5n1}"
puts "Максимальная общая подпоследовательность между COVID-19 (Ухань) и Delta: #{max_common_covid_delta}"
puts "Отношение длины максимальной общей подпоследовательности к общей длине генома: #{ratio_lcs_covid_delta}"
attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
 gttctctaaa cgaactttaa aatctgtgtg gctgtcactc ggctgcatgc ttagtgcact
 cacgcagtat aattaataac taattactgt cgttgacagg acacgagtaa ctcgtctatc
 ttctgcaggc tgcttacggt ttcgtccgtg ttgcagccga tcatcagcac atctaggttt
 cgtccgggtg tgaccgaaag gtaagatgga gagccttgtc cctggtttca acgagaaaac
 acacgtccaa ctcagtttgc ctgttttaca ggttcgcgac gtgctcgtac gtggctttgg
 agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
 cttagtagaa gttgaaaaag gcgttttgcc tcaacttgaa cagccctatg tgttcatcaa
 acgttcggat gctcgaactg cacctcatgg tcatgttatg gttgagctgg tagcagaact
 cgaaggcatt cagtacggtc gtagtggtga gacacttggt gtccttgtcc ctcatgtggg
 cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
 tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
 tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg ttacccgtga
 actcatgcgt gagcttaacg gaggggcata cactcgctat gtcgataaca acttctgtgg
 ccctgatggc taccctcttg agtgcattaa agaccttcta gcacgtgctg gtaaagcttc
 atgcactttg tccgaacaac tggactttat tgacactaag aggggtgtat actgctgccg
 tgaacatgag catgaaattg cttggtacac ggaacgttct gaaaagagct atgaattgca
 gacacctttt gaaattaaat tggcaaagaa atttgacacc ttcaatgggg aatgtccaaa
 ttttgtattt cccttaaatt ccataatcaa gactattcaa ccaagggttg aaaagaaaaa
 gcttgatggc tttatgggta gaattcgatc tgtctatcca gttgcgtcac caaatgaatg
 caaccaaatg tgcctttcaa ctctcatgaa gtgtgatcat tgtggtgaaa cttcatggca
 gacgggcgat tttgttaaag ccacttgcga attttgtggc actgagaatt tgactaaaga
 aggtgccact acttgtggtt acttacccca aaatgctgtt gttaaaattt attgtccagc
 atgtcacaat tcagaagtag gacctgagca tagtcttgcc gaataccata atgaatctgg
 cttgaaaacc attcttcgta agggtggtcg cactattgcc tttggaggct gtgtgttctc
 ttatgttggt tgccataaca agtgtgccta ttgggttcca cgtgctagcg ctaacatagg
 ttgtaaccat acaggtgttg ttggagaagg ttccgaaggt cttaatgaca accttcttga
 aatactccaa aaagagaaag tcaacatcaa tattgttggt gactttaaac ttaatgaaga
 gatcgccatt attttggcat ctttttctgc ttccacaagt gcttttgtgg aaactgtgaa
 aggtttggat tataaagcat tcaaacaaat tgttgaatcc tgtggtaatt ttaaagttac
 aaaaggaaaa gctaaaaaag gtgcctggaa tattggtgaa cagaaatcaa tactgagtcc
 tctttatgca tttgcatcag aggctgctcg tgttgtacga tcaattttct cccgcactct
 tgaaactgct caaaattctg tgcgtgtttt acagaaggcc gctataacaa tactagatgg
 aatttcacag tattcactga gactcattga tgctatgatg ttcacatctg atttggctac
 taacaatcta gttgtaatgg cctacattac aggtggtgtt gttcagttga cttcgcagtg
 gctaactaac atctttggca ctgtttatga aaaactcaaa cccgtccttg attggcttga
 agagaagttt aaggaaggtg tagagtttct tagagacggt tgggaaattg ttaaatttat
 ctcaacctgt gcttgtgaaa ttgtcggtgg acaaattgtc acctgtgcaa aggaaattaa
 ggagagtgtt cagacattct ttaagcttgt aaataaattt ttggctttgt gtgctgactc
 tatcattatt ggtggagcta aacttaaagc cttgaattta ggtgaaacat ttgtcacgca
 ctcaaaggga ttgtacagaa agtgtgttaa gttaaatcca gagaagaaac tggcctactc
 atgcc
 tctaaaagcc ccaaaagaaa ttatcttctt agagggagaa acacttccca cagaagtgtt
 aacagaggaa gttgtcttga aaactggtga tttacaacca ttagaacaac ctactagtga
 agctgttgaa gctccattgg ttggtacacc agtttgtatt aacgggctta tgttgctcga
 aatcaaagac acagaaaagt actgtgccct tgcacctaat atgatggtaa caaacaatac
 cttcacactc aaaggcggtg caccaacaaa ggttactttt ggtgatgaca ctgtgataga
 agtgcaaggt tacaagagtg tgaatatcac ttttgaactt gatgaaagga ttgataaagt
 acttaatgag aagtgctctg cctatacagt tgaactcggt acagaagtaa atgagttcgc
 ctgtgttgtg gcagatgctg tcataaaaac tttgcaacca gtatctgaat tacttacacc
 actgggcatt gatttagatg agtggagtat ggctacatac tacttatttg atgagtctgg
 tgagtttaaa ttggcttcac atatgtattg ttctttctac cctccagatg aggatgaaga
 agaaggtgat tgtgaagaag aagagtttga gccatcaact caatatgagt atggtactga
 agatgattac caaggtaaac ctttggaatt tggtgccact tctgctgctc ttcaacctga
 agaagagcaa gaagaagatt ggttagatga tgatagtcaa caaactgttg gtcaacaaga
 cggcagtgag gacaatcaga caactactat tcaaacaatt gttgaggttc aacctcaatt
 agagatggaa cttacaccag ttgttcagac tattgaagtg aatagtttta gtggttattt
 aaaacttact gacaatgtat acattaaaaa tgcagacatt gtggaagaag ctaaaaaggt
 aaaaccaaca gtggttgtta atgcagccaa tgtttacctt aaacatggag gaggtgttgc
 aggagcctta aataaggcta ctaacaatgc catgcaagtt gaatctgatg attacatagc
 tactaatgga ccacttaaag tgggtggtag ttgtgtttta agcggacaca atcttgctaa
 acactgtctt catgttgtcg gcccaaatgt taacaaaggt gaagacattc aacttcttaa
 gagtgcttat gaaaatttta atcagcacga agttctactt gcaccattat tatcagctgg
 tatttttggt gctgacccta tacattcttt aagagtttgt gtagatactg ttcgcacaaa
 tgtctactta gctgtctttg ataaaaatct ctatgacaaa cttgtttcaa gctttttgga
 aatgaagagt gaaaagcaag ttgaacaaaa gatcgctgag attcctaaag aggaagttaa
 gccatttata actgaaagta aaccttcagt tgaacagaga aaacaagatg ataagaaaat
 caaagcttgt gttgaagaag ttacaacaac tctggaagaa actaagttcc tcacagaaaa
 cttgttactt tatattgaca ttaatggcaa tcttcatcca gattctgcca ctcttgttag
 tgacattgac atcactttct taaagaaaga tgctccatat atagtgggtg atgttgttca
 agagggtgtt ttaactgctg tggttatacc tactaaaaag gctggtggca ctactgaaat
 gctagcgaaa gctttgagaa aagtgccaac agacaattat ataaccactt acccgggtca
 gggtttaaat ggttacactg tagaggaggc aaagacagtg cttaaaaagt gtaaaagtgc
 cttttacatt ctaccatcta ttatctctaa tgagaagcaa gaaattcttg gaactgtttc
 ttggaatttg cgagaaatgc ttgcacatgc agaagaaaca cgcaaattaa tgcctgtctg
 tgtggaaact aaagccatag tttcaactat acagcgtaaa tataagggta ttaaaataca
 agagggtgtg gttgattatg gtgctagatt ttacttttac accagtaaaa caactgtagc
 gtcacttatc aacacactta acgatctaaa tgaaactctt gttacaatgc cacttggcta
 tgtaacacat ggcttaaatt tggaagaagc tgctcggtat atgagatctc tcaaagtgcc
 agctacagtt tctgtttctt cacctgatgc tgttacagcg tataatggtt atcttacttc
 ttcttctaaa acacctgaag aacattttat tgaaaccatc tcacttgctg gttcctataa
 agattggtcc tattctggac aatctacaca actaggtata gaatttctta agagaggtga
 taaaagtgta tattacacta gtaatcctac cacattccac ctagatggtg aagttatcac
 ctttgacaat cttaagacac ttctttcttt gagagaagtg aggactatta aggtgtttac
 aacagtagac aacattaacc tccacacgca agttgtggac atgtcaatga catatggaca
 acagtttggt ccaacttatt tggatggagc tgatgttact aaaataaaac ctcataattc
HEREDOC).slice(0, 6000)

h5n1_sequence = clean_sequence(<<~HEREDOC
agcaaaagca gggtagataa tcactcaatg agtgacatcg aagccatggc atctcaaggc
 accaaacgat catatgaaca aatggagact ggtggggaac gccaggatgc cacagaaatc
 agagcatctg tcggaaaaat gattggtgga atcgggaaat tctacatcca aatgtgcact
 gaactcaaac tcagtgacta tgagggacga ctaatccaaa atagcataac aatagagaga
 atggtgctct ctgcttttga tgagagaaga aataaatacc tagaagagca tcccagtgct
 gggaaggatc ctaagaaaac tggaggaccc atatatagaa gagtagacgg aaagtggatg
 agagaactca ttctttatga caaagaagaa ataaggagag tttggcgcca agcaaataat
 ggtgaagatg caacagctgg tcttactcat atcatgattt ggcattccaa tctgaatgat
 gccacgtacc agagaacaag agcgcttgtt cgcaccggaa tggatcccag aatgtgctct
 ctaatgcaag gttcaacact tcccagaagg tctggggccg caggtgctgc agtgaaagga
 gttggaacaa tagcaatgga attgatcaga atgatcaaac gtgggatcaa tgaccgaaac
 ttctggagag gtgaaaatgg acgaaggaca agggttgcat atgaaagaat gtgcaatatt
 ctcaaaggaa aatttcagac agctgcccag agggcaatga tggatcaagt gagagaaagt
 cggaacccag gaaacgctga gattgaagat ctcattttcc tggcacggtc agcacttatt
 ctaaggggat cagttgcaca taagtcttgc ctgcctgctt gtgtgtatgg gcttgcagtg
 gcaagtgggc atgactttga aagggaaggg tattcactgg tcgggataga cccgtttaaa
 ttactccaaa acagtcaagt gttcagcctg ataagaccaa atgaaaaccc agctcacaag
 agtcaattag tgtggatggc atgccactct gctgcatttg aggatctaag ggtatcaagt
 ttcataagag ggaagaaagt gattccaaga ggaaagcttt ccacaagagg agttcagatt
 gcttcaaatg agaatgtgga agccatggat tccaatacct tagagctgag aagcagatac
 tgggccataa ggaccagaag tggaggaaat accaatcaac aaaaagcatc cgctggccag
 atcagtgtgc aacctacatt ctcagtgcaa cggaatctcc cttttgaaag agcaaccgtt
 atggcagcct tcagcgggaa caatgaagga cggacatccg atatgcgaac agaagttata
 aagatgatgg aaagtgcaaa gccagaagat ttgtccttcc aggggcgggg agtcttcgag
 ctctcggacg aaaaggcaac gagcccgatc gtgccttcct ttgacatgag taatgaagga
 tcttatttct tcggagacaa tgcagaggag tatgacagtt gaggaaaaat acccttgttt
 ctact
HEREDOC).slice(0, 6000)
delta_sequence = clean_sequence(<<~HEREDOC
caaaccaacc aactttcgat ctcttgtaga tctgttctct aaacgaactt taaaatctgt
 gtggctgtca ctcggctgca tgcttagtgc actcacgcag tataattaat aactaattac
 tgtcgttgac aggacacgag taactcgtct atcttctgca ggctgcttac ggtttcgtcc
 gttttgcagc cgatcatcag cacatctagg ttttgtccgg gtgtgaccga aaggtaagat
 ggagagcctt gtccctggtt tcaacgagaa aacacacgtc caactcagtt tgcctgtttt
 acaggttcgc gacgtgctcg tacgtggctt tggagactcc gtggaggagg tcttatcaga
 ggcacgtcaa catcttaaag atggcacttg tggcttagta gaagttgaaa aaggcgtttt
 gcctcaactt gaacagccct atgtgttcat caaacgttcg gatgctcgaa ctgcacctca
 tggtcatgtt atggttgagc tggtagcaga actcgaaggc attcagtacg gtcgtagtgg
 tgagacactt ggtgtccttg tccctcatgt gggcgaaata ccagtggctt accgcaaggt
 tcttcttcgt aagaacggta ataaaggagc tggtggccat agttacggcg ccgatctaaa
 gtcatttgac ttaggcgacg agcttggcac tgatccttat gaagattttc aagaaaactg
 gaacactaaa catagcagtg gtgttacccg tgaactcatg cgtgagctta acggaggggc
 atacactcgc tatgtcgata acaacttctg tggccctgat ggctaccctc ttgagtgcat
 taaagacctt ctagcacgtg ctggtaaagc ttcatgcact ttgtccgaac aactggactt
 tattgacact aagaggggtg tatactgctg ccgtgaacat gagcatgaaa ttgcttggta
 cacggaacgt tctgaaaaga gctatgaatt gcagacacct tttgaaatta aattggcaaa
 gaaatttgac accttcaatg gggaatgtcc aaattttgta tttcccttaa attccataat
 caagactatt caaccaaggg ttgaaaagaa aaagcttgat ggctttatgg gtagaattcg
 atctgtctat ccagttgcgt caccaaatga atgcaaccaa atgtgccttt caactctcat
 gaagtgtgat cattgtggtg aaacttcatg gcagacgggc gattttgtta aagccacttg
 cgaattttgt ggcactgaga atttgactaa agaaggtgcc actacttgtg gttacttacc
 ccaaaatgct gttgttaaaa tttattgtcc agcatgtcac aattcagaag taggacctga
 gcatagtctt gccgaatacc ataatgaatc tggcttgaaa accattcttc gtaagggtgg
 tcgcactatt gcctttggag gctgtgtgtt ctcttatgtt ggttgccata acaagtgtgc
 ctattgggtt ccacgtgcta gcgctaacat aggttgtaac catacaggtg ttgttggaga
 aggttccgaa ggtcttaatg acaaccttct tgaaatactc caaaaagaga aagtcaacat
 caatattgtt ggtgacttta aacttaatga agagatcgcc attattttgg catctttttc
 tgcttccaca agtgcttttg tggaaactgt gaaaggtttg gattataaag cattcaaaca
 aattgttgaa tcctgtggta attttaaagt tacaaaagga aaagctaaaa aaggtgcctg
 gaatattggt gaacagaaat caatactgag tcctctttat gcatttgcat cagaggctgc
 tcgtgttgta cgatcaattt tctcccgcac tcttgaaact gctcaaaatt ctgtgcgtgt
 tttacagaag gccgctataa caatactaga tggaatttca cagtattcac tgagactcat
 tgatgctatg atgttcacat ctgatttggc tactaacaat ctagttgtaa tggcctacat
 tacaggtggt gttgttcagt tgacttcgca gtggctaact aacatctttg gcactgttta
 tgaaaaactc aaacccgtcc ttgattggct tgaagagaag tttaaggaag gtgtagagtt
 tcttagagac ggttgggaaa ttgttaaatt tatctcaacc tgtgcttgtg aaattgtcgg
 tggacaaatt gtcacctgtg caaaggaaat taaggagagt gttcagacat tctttaagct
 tgtaaataaa tttttggctt tgtgtgctga ctctatcatt attggtggag ctaaacttaa
 agccttgaat ttaggtgaaa catttgtcac gcactcaaag ggattgtaca gaaagtgtgt
 taaatccaga gaagaaactg gcctactcat gcctctaaaa gccccaaaag aaattatctt
 cttagaggga gaaacacttc ccacagaagt gttaacagag gaagttgtct tgaaaactgg
 tgatttacaa ccattagaac aacctactag tgaagctgtt gaagctccat tggttggtac
 accagtttgt attaacgggc ttatgttgct cgaaatcaaa gacacagaaa agtactgtgc
 ccttgcacct aatatgatgg taacaaacaa taccttcaca ctcaaaggcg gtgcaccaac
 aaaggttact tttggtgatg acactgtgat agaagtgcaa ggttacaaga gtgtgaatat
 cacttttgaa cttgatgaaa ggattgataa agtacttaat gagaagtgct ctgcctatac
 agttgaactc ggtacagaag taaatgagtt cgcctgtgtt gtggcagatg ctgtcataaa
 aactttgcaa ccagtatctg aattacttac accactgggc attgatttag atgagtggag
 tatggctaca tactacttat ttgatgagtc tggtgagttt aaattggctt cacatatgta
 ttgttctttt taccctccag atgaggatga agaagaaggt gattgtgaag aagaagagtt
 tgagccatca actcaatatg agtatggtac tgaagatgat taccaaggta aacctttgga
 atttggtgcc acttctgctg ctcttcaacc tgaagaagag caagaagaag attggttaga
 tgatgatagt caacaaactg ttggtcaaca agacggcagt gaggacaatc agacaactac
 tattcaaaca attgttgagg ttcaacctca attagagatg gaacttacac cagttgttca
 gactattgaa gtgaatagtt ttagtggtta tttaaaactt actgacaatg tatacattaa
 aaatgcagac attgtggaag aagctaaaaa ggtaaaacca acagtggttg ttaatgcagc
 caatgtttac cttaaacatg gaggaggtgt tgcaggagcc ttaaataagg ctactaacaa
 tgccatgcaa gttgaatctg atgattacat agctactaat ggaccactta aagtgggtgg
 tagttgtgtt ttaagcggac acaatcttgc taaacactgt cttcatgttg tcggcccaaa
 tgttaacaaa ggtgaagaca ttcaacttct taagagtgct tatgaaaatt ttaatcagca
 cgaagttcta cttgcaccat tattatcagc tggtattttt ggtgctgacc ctatacattc
 tttaagagtt tgtgtagata ctgttcgcac aaatgtctac ttagctgtct ttgataaaaa
 tctctatgac aaacttgttt caagcttttt ggaaatgaag agtgaaaagc aagttgaaca
 aaagatcgct gagattccta aagaggaagt taagccattt ataactgaaa gtaaaccttc
 agttgaacag agaaaacaag atgataagaa aatcaaagct tgtgttgaag aagttacaac
 aactctggaa gaaactaagt tcctcacaga aaacttgtta ctttatattg acattaatgg
 caatcttcat ccagattctg ccactcttgt tagtgacatt gacatcactt tcttaaagaa
 agatgctcca tatatagtgg gtgatgttgt tcaagagggt gttttaactg ctgtggttat
 acctactaaa aagtctggtg gcactactga aatgctagcg aaagctttga gaaaagtgcc
 aacagacaat tatataacca cttacccggg tcagggttta aatggttaca ctgtagagga
 ggcaaagaca gtgcttaaaa agtgtaaaag tgccttttac attctaccat ctattatctc
 taatgagaag caagaaattc ttggaactgt ttcttggaat ttgcgagaaa tgcttgcaca
 tgcagaagaa acacgcaaat taatgcctgt ctgtgtggaa actaaagcca tagtttcaac
 tatacagcgt aaatataagg gtattaaaat acaagagggt gtggttgatt atggtgctag
 attttacttt tacaccagta aaacaactgt agcgtcactt atcaacacac ttaacgatct
 aaatgaaact cttgttacaa tgccacttgg ctatgtaaca catggcttaa atttggaaga
 agctgctcgg tatatgagat ctctcaaagt gccagctaca gtttctgttt cttcacctga
 tgctgttaca gcgtataatg gttatcttac ttcttcttct aaaacacctg aagaacattt
 tattgaaacc atctcacttg ctggttccta taaagattgg tcctattctg gacaatctac
 acaactaggt atagaatttc ttaagagagg tgataaaagt gtatattaca ctagtaatcc
 taccacattc cacctagatg gtgaagttat cacctttgac aatcttaaga cacttctttc
 tttgagagaa gtgaggacta ttaaggtgtt tacaacagta gacaacatta acctccacac
 gcaagttgtg gacatgtcaa tgacatatgg acaacagttt ggtccaactt atttggatgg
 agctgatgtt actaaaataa aacctcataa ttcacatgaa ggtaaaacat tttatgtttt
HEREDOC).slice(0, 6000)

covid_cleaned = clean_sequence(covid_sequence)
h5n1_cleaned = clean_sequence(h5n1_sequence)
delta_cleaned = clean_sequence(delta_sequence)

min_unique_covid_h5n1 = find_min_unique(covid_cleaned, h5n1_cleaned)
min_unique_h5n1_covid = find_min_unique(h5n1_cleaned, covid_cleaned)

min_unique_covid_h5n1_delta = find_min_unique(covid_cleaned, h5n1_cleaned + delta_cleaned)
min_unique_delta_covid_h5n1 = find_min_unique(delta_cleaned, h5n1_cleaned + covid_cleaned)

max_common_covid_delta, ratio_lcs_covid_delta = find_lcs(covid_cleaned, delta_cleaned)

puts "Минимальная специфичная последовательность для COVID-19 (Ухань), отсутствующая в H5N1: #{min_unique_covid_h5n1}"
puts "Минимальная специфичная последовательность для H5N1, отсутствующая в COVID-19 (Ухань): #{min_unique_h5n1_covid}"
puts "Минимальная специфичная последовательность для COVID-19 (Ухань), отсутствующая в Delta и H5N1: #{min_unique_covid_h5n1_delta}"
puts "Минимальная специфичная последовательность для Delta, отсутствующая в COVID-19 (Ухань) и H5N1: #{min_unique_delta_covid_h5n1}"
puts "Максимальная общая подпоследовательность между COVID-19 (Ухань) и Delta: #{max_common_covid_delta}"
puts "Отношение длины максимальной общей подпоследовательности к общей длине генома: #{ratio_lcs_covid_delta}"
