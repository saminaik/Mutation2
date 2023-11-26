package com.example.mutation2;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static java.lang.Math.max;
import static org.junit.jupiter.api.Assertions.*;

class Mutation2ApplicationTest {

    @Test
    void main() {
    }

    @Test
    void add() {
        Assertions.assertEquals(new Mutation2Application().add(5,2),7);
    }
    @Test
    void sub(){
        Assertions.assertEquals(new Mutation2Application().sub(10,5),5);
    }
    Mutation2Application obj=new Mutation2Application();
    @Test
    public void TestEditDistance() {
        String str1 = "sunday";
        String str2 = "saturday";

        int n = str1.length(), m = str2.length();
        int[][] dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);

        assertEquals(3, obj.EditDistance(str1, str2, n, m, dp), 0.0);

        str1 = "dinitrophenylhydrazine";
        str2 = "benzalphenylhydrazone";

        n = str1.length();
        m = str2.length();
        dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);
        assertEquals(7, obj.EditDistance(str1, str2, n, m, dp), 0.0);

        str1 = "quantumcomputingdynamicprogramming";
        str2 = "quancmpdmicprmmi";

        n = str1.length();
        m = str2.length();
        dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);
        assertEquals(18, obj.EditDistance(str1, str2, n, m, dp), 0.0);

    }

    @Test
    public void TestKMPAlgorithm()
    {
        String txt = "ABABDABACDABABCABAB";
        String pat = "ABABCABAB";
        assertEquals(10, obj.KMPSearch(pat, txt));

        txt = "HSKCMJEHSTABRABCUHSGDFYUGSUDVBDUSHVBDUHBVSUHDG";
        pat = "HSKC";
        assertEquals(0, obj.KMPSearch(pat, txt));

        txt = "BDUHBVSUHDG";
        pat = "DJKNSJK";
        assertEquals(-1, obj.KMPSearch(pat, txt));

        txt = "ASDFGH";
        pat = "GHI";
        assertEquals(-1, obj.KMPSearch(pat, txt));
    }

    @Test
    public void TestRabinKarp()
    {
        String txt = "HSKCMJEHSTABRABCUHSGDFYUGSUDVBDUSHVBDUHBVSUHDG";
        String pat = "ABRA";
        // A prime number
        int q = 101;
        assertEquals(10, obj.RabinKarp(pat, txt, q), 0.0);

        txt = "JHDYEKCHSPLQJAYDABRANCJKSDNFNJKDNCCDJKNFJKDN";
        pat = "JHDY";
        assertEquals(0, obj.RabinKarp(pat, txt, q), 0.0);

        txt = "BDUHBVSUHDG";
        pat = "DJKNSJK";
        assertEquals(-1, obj.RabinKarp(pat, txt, q), 0.0);

        txt = "BDUHBVSUHDG";
        pat = "DGFDHBJDKBJ";
        assertEquals(-1, obj.RabinKarp(pat, txt, q), 0.0);

        txt = "dhjbdjs";
        pat = "js";
        assertEquals(5, obj.RabinKarp(pat, txt, q), 0.0);
    }

    @Test
    public void TestZAlgorithm()
    {
        String text = "JHDYEKCHSPLQJAYDABRANCJKSDNFNJKDNCCDJKNFJKDN";
        String pattern = "DABRA";
        assertEquals(15, obj.ZAlgorithm(text, pattern), 0.0);

        text = "HSKCMJEHSTABRABCUHSGDFYUGSUDVBDUSHVBDUHBVSUHDG";
        pattern = "ABRA";
        assertEquals(10, obj.ZAlgorithm(text, pattern), 0.0);

        text = "BDUHBVSUHDG";
        pattern = "DJKNSJKJSDJWESUJSD";
        assertEquals(-1, obj.ZAlgorithm(text, pattern), 0.0);

        text = "";
        pattern = "dfsdfs";
        assertEquals(-1, obj.ZAlgorithm(text, pattern), 0.0);

        text = "ABCDFFGH";
        pattern = "H";
        assertEquals(7, obj.ZAlgorithm(text, pattern), 0.0);

        text = "aaaaaaaaaaaa";
        pattern = "aa";
        assertEquals(0, obj.ZAlgorithm(text, pattern));
    }

    @Test
    public void TestShortestCommonSequence()
    {
        String X = "AGGTB";
        String Y = "GXTXAYB";

        int[][] lookup
                = new int[X.length() + 1][Y.length() + 1];

        assertEquals(9, obj.superSeq(X, Y, X.length(), Y.length(), lookup), 0.0);

        X = "apqrstu";
        Y = "kplrmntuo";
        lookup
                = new int[X.length() + 1][Y.length() + 1];
        assertEquals(12, obj.superSeq(X, Y, X.length(), Y.length(), lookup), 0.0);
    }

    @Test
    public void TestLCS()
    {
        String s1 = "AGGTAB";
        String s2 = "GXTXAYB";

        char[] X=s1.toCharArray();
        char[] Y=s2.toCharArray();
        int m = X.length;
        int n = Y.length;
        assertEquals(4, obj.LCS(X, Y, m, n), 0.0);

        s1 = "ABRACADABRA";
        s2 = "YABBADABBADOO";
        X=s1.toCharArray();
        Y=s2.toCharArray();
        m = X.length;
        n = Y.length;
        assertEquals(7, obj.LCS(X, Y, m, n), 0.0);
    }

    @Test
    public void TestLSFP()
    {
        String s = "adbabd";

        // Length of given string
        int n = s.length();

        // Function call
        assertEquals(6, obj.longestSubstring(s, n), 0.0);
        s = "pqrstuabcdefffghijkahgjgujgjudykisykiuljiktytertersrgdghfhfhfhxf";
        n = s.length();
        assertEquals(9, obj.longestSubstring(s, n), 0.0);
    }

    @Test
    public void TestLVP()
    {
        assertEquals(8, obj.LVP("((()()()()(((())", 16), 0.0);
        assertEquals(0, obj.LVP("((((((", 6), 0.0);
        assertEquals(2, obj.LVP("())))))))", 9), 0.0);
        assertEquals(0, obj.LVP(")", 1));
        assertEquals(0, obj.LVP("", 0));
    }

    @Test
    public void TestLongestCommonPrefix() {
        String[] input = {"geeksforgeeks", "geeks", "geek", "geezer"};
        assertEquals("gee", obj.longestCommonPrefix(input));

        String[] input1 = {"zacchsdjiwsbd", "zaccvbhdbfeui", "zaccbuivghdreuhg", "zaccfguieh", "zaccbfuid", "zaccghuvbifryhe", "zaccvbuidfb", "zacchbfui", "zaccbuidf"};
        assertEquals("zacc", obj.longestCommonPrefix(input1));

        String[] input2 = {"reflower","flow","flight"};
        assertEquals("", obj.longestCommonPrefix(input2));

        String[] input3 = {"reflower"};
        assertEquals("reflower", obj.longestCommonPrefix(input3));

        String[] input4 = {};
        assertEquals("", obj.longestCommonPrefix(input4));

        String[] input5 = {"dcnsj", "dsdsd", "djsjs"};
        assertEquals("d", obj.longestCommonPrefix(input5));
    }

    @Test
    public void TestLongestPalindromicSubsequence() {
        String seq = "GEEKSFORGEEKS";
        int n = seq.length();
        assertEquals(5, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);

        String text = "abbcbbaj";
        n = text.length();
        assertEquals(7, obj.lps(text.toCharArray(), 0, n - 1), 0.0);

//        seq = "bdusiwbauibwuisdbcwquibduiwboqibxwuqibxquiwbx";
//        n = seq.length();
//        assertEquals(23, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);
    }

    @Test
    public void TestManachersAlgorithm() {
        String text = "abacdfgdcaba";
        assertEquals("aba", obj.ManachersAlgorithm(text));
        text = "abcdeghghijkalmnahgyuwdggggaygufgsiugfisuuuuaaagggauuuuuuuhhhabbbbiiiiiiiiisssssssssoooooooaaauuuuuu";
        assertEquals("iiiiiiiii", obj.ManachersAlgorithm(text));

        text = "jdfhsjkhdf";
        assertEquals("d", obj.ManachersAlgorithm(text));

        text = "aabbbbaababababa";
        assertEquals("ababababa", obj.ManachersAlgorithm(text));
    }

    @Test
    public void TestBoyerMoore() {
        char txt[] = "ABAAABCD".toCharArray();
        char pat[] = "ABC".toCharArray();
        assertEquals(4, obj.BoyerMoore(txt, pat));

        txt = "BDUHBVSUHDG".toCharArray();
        pat = "DJKNSJKJSDJWESUJSD".toCharArray();
        assertEquals(-1, obj.BoyerMoore(txt, pat), 0.0);
    }

    @Test
    public void TestSequenceAlignment() {
        String gene1 = "AGGGCT";
        String gene2 = "AGGCA";

        // initialising penalties
        // of different types
        int misMatchPenalty = 3;
        int gapPenalty = 2;

        // calling the function to
        // calculate the result
        assertEquals(5, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "CG";
        gene2 = "CA";
        misMatchPenalty = 3;
        gapPenalty = 7;
        assertEquals(3, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "";
        gene2 = "AGGCA";
        misMatchPenalty = 3;
        gapPenalty = 2;
        assertEquals(10, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "AGGGCT";
        gene2 = "";
        misMatchPenalty = 3;
        gapPenalty = 2;
        assertEquals(12, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "sdbsjkdnbjksbndjk";
        gene2 = "snajksnjkansjkasn";
        misMatchPenalty = 3;
        gapPenalty = 2;
        assertEquals(27, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));


    }

    @Test
    public void TestWildcardPattern() {
        String str = "baaabab";
        String pattern = "*****ba*****ab";
        assertEquals(true, obj.WildcardPattern(str, pattern, str.length(), pattern.length()));

        str = "babaaababaabababbbbbbaabaabbabababbaababbaaabbbaaab";
        pattern = "***bba**a*bbba**aab**";
        assertEquals(true, obj.WildcardPattern(str, pattern, str.length(), pattern.length()));

        str = "";
        pattern = "dksbndfjk";
        assertEquals(false, obj.WildcardPattern(str, pattern, str.length(), pattern.length()));

        str="cojkdsmcos";
        pattern = "";
        assertEquals(false, obj.WildcardPattern(str, pattern, str.length(), pattern.length()));

        str = "";
        pattern = "";
        assertEquals(true, obj.WildcardPattern(str, pattern, str.length(), pattern.length()));
    }

    @Test
    public void TestminPalPartition() {
        String str = "ababbbabbababa";
        assertEquals(3, obj.minPalPartition(str));
    }

    @Test
    public void TestSparseSearch() {
        String arr[] = {"for", "hjcbs", "", "", "", "", "ide",
                "pracs", "", "", "", "quiz"};
        String x = "for";
        int n = x.length();
        assertEquals(0, obj.SparseSearch(arr, x, n));

        String arr1[] = {"amdck", "bdhsd", "", "", "", "", "ide"};
        String x1 = "bdhsd";
        int n1 = x1.length();
        assertEquals(1, obj.SparseSearch(arr1, x1, n1));


        String arr2[] = {"for", "hjcbs", "", "", "", "cnjdk"};
        String x2 = "jki";
        int n2 = x2.length();
        assertEquals(-1, obj.SparseSearch(arr2, x2, n2));

        String arr3[] = {"", "", ""};
        String x3 = "jki";
        int n3 = x3.length();
        assertEquals(-1, obj.SparseSearch(arr3, x3, n3));

        String arr4[] = {"", "", "", "dwe", "", "zz"};
        String x4 = "zz";
        int n4 = x4.length();
        assertEquals(-1, obj.SparseSearch(arr4, x4, n4));
    }

    @Test
    public void TestLongestRepeatingSubSeq() {
        String str = "aabb";
        assertEquals(2, obj.LongestRepeatingSubSeq(str));
    }

    @Test
    public void TestLongestPrefixSuffix() {
        String str = "abcab";
        assertEquals(2, obj.longestPrefixSuffix(str));
        str = "asdfgjndjkdnjkendjkvnsdjknjkdnwuhnclwiosdcnklwnasdfg";
        assertEquals(5, obj.longestPrefixSuffix(str));
        str = "fcdjksncjksndcjksd";
        assertEquals(0, obj.longestPrefixSuffix(str));
    }

    @Test
    public void TestKVowelWords() {
        int N = 3;
        int K = 3;
        assertEquals(17576, obj.KVowelWords(N, K));
    }

    @Test
    public void TestLRR() {
        String str = "ncbjknsdjkcnsjkancjksdncjksdncjk";
        assertEquals("jknsdjkcnsjkancjksdncjksdncjkncb", obj.leftrotate(str, 3));
        assertEquals("jkncbjknsdjkcnsjkancjksdncjksdnc", obj.rightrotate(str, 2));

    }

    @Test
    public void TestReverseVowel() {
        String str= "hello world";
        assertEquals("hollo werld", obj.reverseVowel(str));

    }

    @Test
    public void TestRepeatedStringMatch() {
        String a = "abcd", b = "cdabcdab";
        assertEquals(3, obj.repeatedStringMatch(a, b));
        a = "abcde";
        b = "bcdeabcdeabcdeabcdeabcd";
        assertEquals(5, obj.repeatedStringMatch(a, b));
        a = "cdmnslkd";
        b = "iwomeiofnmwo";
        assertEquals(-1, obj.repeatedStringMatch(a, b));
    }


}