//
//  subcipher.h
//  LiveWrite
//
//  Created by Ben Mildenhall on 11/18/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#ifndef LiveWrite_subcipher_h
#define LiveWrite_subcipher_h

/*
 * File: spellcheck.c
 * Author: YOUR NAME HERE
 * ----------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <time.h>
#include <math.h>
using namespace std;

#define CORPUS_SIZE_GUESS 10000
#define MISSPELLINGS_SIZE_GUESS 50
#define MAX_WORD_SIZE 25
#define SUGGESTION_NUM 5

const float kNExistPenalty = 0.f;
const float kWordBonus = 1.5f;


// COPYRIGHT CS 107
bool ReadOneWord(FILE *fp, char buf[])
{
    while (true) { // keep reading until good word or EOF
        if (fscanf(fp, " %25[a-zA-Z]", buf) == 1) { // read at most max letters
            int ch = getc(fp);                      // peek at next char
            if (isspace(ch) || ch == EOF) return true; // accept if end space/EOF
            ungetc(ch, fp);                         // otherwise, put char back
        }
        if (fscanf(fp,"%*s") == EOF) return false; //%*s reads rest of token and discards
    }
}

void ToLowerCase(char *str) {
    int length = strlen(str);
    //To ensure lower case, flip on the bit in the 2^5 slot, which
    //controls lowercase-ness
    for (int i = 0; i < length; i++) str[i] |= 0x20;
}

void FillInCorpus(FILE *fp, map<string, int> &freqs) {
    char buf[26];
    while (ReadOneWord(fp, buf)) {
    	ToLowerCase(buf);
        freqs[string(buf)]++;
    }
}

int Min3i(int a, int b, int c) {
    return min(a, min(b, c));
}

int EditDistance(const char *word, const char *dest, int dist, int max) {
    if (dist > max) return dist;
    if (word[0] == dest[0]) return word[0] ? EditDistance(word + 1, dest + 1, dist, max) : dist;
    if (word[0] == '\0') return EditDistance(word, dest + 1, dist + 1, max);
    if (dest[0] == '\0') return EditDistance(word + 1, dest, dist + 1, max);
    return Min3i(EditDistance(word, dest + 1, dist + 1, max),
                 EditDistance(word + 1, dest, dist + 1, max),
                 EditDistance(word + 1, dest + 1, dist + 1, max));
}

string SpellCheck(map<string, int> freqs, string word) {
    int maxdist = INT_MAX;
    int bestfreq;
    string bestSug;
    //cout << "Spell checking " << word << "..." << endl;
    for (auto &p : freqs) {
        int dist = EditDistance(word.c_str(), p.first.c_str(), 0, maxdist);
        if (dist < maxdist || (dist == maxdist && p.second < bestfreq)) {
            maxdist = dist;
            bestfreq = p.second;
            bestSug = p.first;
            //cout << "New best is " << bestSug << " with edit dist " << maxdist << endl;
        }
    }
    return bestSug;
}

void RandomPermutation(vector<int> &vals) {
    for (int i = 0; i < vals.size() - 1; ++i) {
        int index = i + (rand() % (vals.size() - i));
        swap(vals[i], vals[index]);
    }
}

void AlphaPermute(vector<int> &vals) {
    vals.resize(26);
    for (int i = 0; i < 26; ++i)
        vals[i]  = i;
    RandomPermutation(vals);
}

string PermuteString(vector<int> &key, string str) {
    string result = str;
    for (char &c : result)
        if (c >= 'a' && c <= 'z')
            c = key[c - 'a'] + 'a';
    return result;
}

void PermuteStringVector(vector<int> &key, vector<string> &strings) {
    for (string &s : strings)
        s = PermuteString(key, s);
}

void GetNGrams(map<string, int> &wordFreqs, map<string, int> &ngrams, int n) {
    for (auto &p : wordFreqs) {
        string word = string('$', n-1) + p.first + string('^', n-1);
        int count = p.second;
        for (int i = 0; i < word.size() - n; ++i)
            ngrams[word.substr(i, n)] += count;
    }
}

void RandomSwap(vector<int> &vals) {
    if (vals.size() <= 1) return;
    int a = rand() % vals.size();
    int b = rand() % vals.size();
    while (a == b)
        b = rand() % vals.size();
    swap(vals[a], vals[b]);
}

float NGramScore(map<string, int> &ngrams, vector<string> &input) {
    int n = ngrams.begin()->first.length();
    float score = 0.f;
    for (string &s : input)
        for (int i = 0; i < s.size() - n; ++i)
            score += log(1.f + ngrams[s.substr(i, n)]);
    return score;
}

float WordScore(map<string, int> &freqs, vector<string> &input, int n) {
    float score = 0.f;
    for (string &s : input) {
        string curr = s.substr(n - 1, s.length() - 2 * (n - 1));
        if (freqs.find(curr) == freqs.end())
            score += kNExistPenalty;
        else
            score += log((float)freqs[curr]);
    }
    return score * kWordBonus;
}

void RunInference(int trials, int iters, map<string, int> &freqs,
                  map<string, int> &ngrams, vector<string> &input, int n) {
    vector<string> modIn(input);
    for (string &s : modIn)
        s = string(n - 1, '$') + s + string(n - 1, '^');
    for (string &s : modIn)
        cout << s << " ";
    float beta = 1.f;
    vector<int> key;
    AlphaPermute(key);
    float bestScore = -INFINITY;
    vector<int> bestKey;
    for (int t = 0; t < trials; ++t) {
        //cout << "Starting trial " << t << endl;
        float oldscore = 0.f;
        AlphaPermute(key);
        for (int i = 0; i < iters; ++i) {
            int swaps = rand() % 3;
            vector<int> newKey(key);
            
            for (int s = 0; s < swaps; ++s)
                RandomSwap(newKey);
            
            swap(newKey[i % 26], newKey[rand() % 26]);
            vector<string> permutation(modIn);
            PermuteStringVector(newKey, permutation);
            float score = NGramScore(ngrams, permutation);
            score += WordScore(freqs, permutation, n);
            if (log((float)rand() / INT_MAX) < beta * (score - oldscore)) {
                oldscore = score;
                key = newKey;
                if (score > bestScore) {
                    bestScore = score;
                    bestKey = key;
                    cout << "New best at " << score << " on iter " << i
                    << ", trial " << t << endl;
                }
            }
        }
    }
    //cout << "Best overall was scored " << bestScore << endl;
    PermuteStringVector(bestKey, input);
    for (string &s : input)
        cout << s << " ";
    PermuteStringVector(bestKey, modIn);
    cout << endl << NGramScore(ngrams, modIn) << " is the ngrams score" << endl;
    cout << WordScore(freqs, modIn, n) << " is the words score" << endl;
    cout << "Returning the above" << endl;
}

void ReadCorpus(string filename, map<string, int> &freqs) {
    FILE *fp = fopen(filename.c_str(), "r");
    if (!fp) {
    	printf("Cannot open file \"%s\"\n", filename.c_str());
        return;
    }
    FillInCorpus(fp, freqs);
    fclose(fp);
    cout << "Successfully read in file " << filename.c_str() << endl;
}



/*
int main(int argc, const char *argv[])
{
    int n = 3;
    srand (time(NULL));
    vector<string> input;
    for (int i = 2; i < argc; ++i)
        input.push_back(string(argv[i]));
    for (string &s : input)
        for (char &c : s)
            if (c >= 'A' && c <= 'Z')
                c += 'z' - 'Z';
    
    
    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
    	printf("Cannot open file \"%s\"\n", argv[1]);
    	return 0;
    }
    map<string, int> freqs;
    FillInCorpus(fp, freqs);
    fclose(fp);
    cout << "Found " << freqs.size() << " words" << endl;
    cout << "tried: " << freqs["tried"] << endl;
    
    int total = 0;
    for (auto &p : freqs) total += p.second;
    cout << "Found " << total << " entries" << endl;
    map<string, int> ngrams;
    GetNGrams(freqs, ngrams, n);
    
    for (string &s : input) cout << s << " ";
    cout << endl;
    
    RunInference(50, 5000, freqs, ngrams, input, n);
    return 0;
}
*/

#endif
