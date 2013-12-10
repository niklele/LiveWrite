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
 #include <fstream>
 #include <iomanip>
 #include <sstream>
using namespace std;

#define CORPUS_SIZE_GUESS 10000
#define MISSPELLINGS_SIZE_GUESS 50
#define MAX_WORD_SIZE 25
#define SUGGESTION_NUM 5

const float kNExistPenalty = 0.f;
const float kWordBonus = 1.5f;
const float epsilon = .00001f;
const static int FREQ_CUTOFF = 0;

static map<string, int> freqs;
static map<string, int> ngrams;
static map<string, vector<string>> cipherCodes;

float randuni() {
    return (float)rand() / INT_MAX;
}

static int trials = 0;
static float biastrialpos = 0.f;
bool trial(float f) {
    bool ret = (randuni() <= f);
    ++trials;
    if (ret) biastrialpos += 1.f / f;
    return ret;
}

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
    string result(str);
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
        string word = string(n-1, '$') + p.first + string(n-1, '$');
        int count = p.second;
        for (int i = 0; i < word.size() - n; ++i) 
            ngrams[word.substr(i, n)] += count;
    }
} 

string CipherCode(const string &s) {
    string code = "";
    int pos = 0, cused = 0;
    for (char c : s) {
        int f = s.find(c);
        if (f == pos++) 
            code += (char)('a' + cused++);
        else 
            code += code[f];
    }
    return code;
}

vector<string> MaskChars(const string &mask, const vector<string>& words) {
    vector<string> retv;
    for (const string &word : words) {
        bool fits = true;
        for (int i = 0; i < mask.length() && fits; ++i)
            fits = (mask[i] == '.' || mask[i] == word[i]);
        if (fits)
            retv.push_back(word);
    }
    return retv;
}

vector<string> MaskCharsWithBanned(const string &mask, const string &banned,
                                                     const vector<string>& words) {
    vector<string> retv;
    for (const string &word : words) {
        bool fits = true;
        for (int i = 0; i < mask.length() && fits; ++i) 
            fits = (mask[i] == word[i] || (mask[i] == '.' 
                                && banned.find(word[i]) == string::npos));
        if (fits)
            retv.push_back(word);
    }
    return retv;
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
        s = string(n - 1, '$') + s + string(n - 1, '$');
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
               
            if (rand() % 2) 
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
        for (int a = 0; a < 25; ++a) {
            for (int b = a + 1; b < 26; ++b) {
                vector<int> newKey(key);
                swap(newKey[a], newKey[b]);
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
                        cout << "New best at " << score << " on special swap " << a << ", " << b
                            << ", trial " << t << endl;
                    }
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

void ReadStringIntMap(const string &filename, map<string, int> &simap, bool ints = true) {
    ifstream infile(filename);
    if (infile.fail()) return;
    while (!infile.eof()) {
        string s;
        int i = 1;
        infile >> s;
        if (ints) infile >> i;
        while ((s[0] > 'z' || s[0] < 'a') && !s.empty()) s = s.substr(1);
        if (s.empty()) continue;
        simap[s] = i;
    }
    infile.close();
}

void GetMatchMat(string &str, vector<float> &mat) {
    int n = str.size();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            mat[j * n + i] = (str[i] == str[j] ? 1.f : 0.f);
        }
    }
}

float CompareMat(vector<float> &m1, vector<float> &m2) {
    int n = min(m1.size(), m2.size());
    float ret = 0;
    for (int i = 0; i < n; ++i) 
        ret += (m1[i] - m2[i]) * (m1[i] - m2[i]) / n / n;
    if (ret < epsilon) ret = epsilon;
    return ret;
}


float NGramScore(vector<string> &input) {
    int n = 3;
    float score = 0.f;
    for (string s : input) {
        s = "$$" + s + "$$";
        for (int i = 0; i < s.size() - n; ++i) 
            score += log(1.f + ngrams[s.substr(i, n)]);
    }
    return score;
}

float WordScore(vector<string> &input) {
    float score = 0.f;
    for (string &s : input) {
        if (freqs.find(s) == freqs.end()) score += kNExistPenalty;
        else score += log((float)freqs[s]);
    }
    return score;
}

void ProbInf(vector<float> &probmat, int nchar, vector<float> &wordL, int threadid = 0) {
    const int iters = 1000000;
    const float beta = .25f;
    const float weights[] = {.001f, 5.f, 1.f};
    vector<float> matchmat(nchar * nchar);
    string guess = "";
    for (int i = 0; i < nchar; ++i) 
        guess += (char)('a' + (rand() % 26));
    
    float score, oldscore = -INFINITY, bestscore = oldscore;
    string bestguess, oldguess = guess;
    vector<string> wordinput;
    int prop, acc;
    prop = acc = 0;
    for (int iter = 0; iter < iters; ++iter) {
        guess = oldguess;
        int numchange = 1 + rand() % 3;
        if (trial(.1f)) numchange = 10;
        for (int z = 0; z < numchange; ++z) {
            int seedind = rand() % nchar;
            int swapind = rand() % nchar;
            char seed = guess[seedind], swap = guess[swapind];
            char ose = seed;
            if (trial(.25f)) {
                seed = 'a' + rand() % 26;
                while (guess.find(seed) != guess.size() && trial(.95f))
                    seed = 'a' + rand() % 26;
                for (int i = 0; i < nchar; ++i) {
                    if (trial(probmat[seedind * nchar + i]))
                        guess[i] = seed;
                }
            } else {
                for (int i = 0; i < nchar; ++i) {
                    if (trial(probmat[seedind * nchar + i]))
                        guess[i] = swap;
                    if (trial(probmat[swapind * nchar + i]))
                        guess[i] = seed;
                }
            }
        }

        GetMatchMat(guess, matchmat);
        float score = weights[0] / CompareMat(probmat, matchmat);
        wordinput.clear();
        int ind = 0;
        for (float l : wordL) {
            wordinput.push_back(guess.substr(ind, l));
            ind += l;
        }

        score += WordScore(wordinput) * weights[1];
        score += NGramScore(wordinput) * weights[2];
        ++prop;
        if (log(randuni()) < (score - oldscore) * beta) {
            ++acc;
        //cout << "OLD " << oldguess << " NEW " << guess << endl;
            oldscore = score;
            oldguess = guess;
            if (bestscore < oldscore) {
                bestscore = oldscore;
                bestguess = oldguess;
                cout << "New best at iter " << iter;
                cout << ", scored " << bestscore << ", guess is: " << endl;
                wordinput.clear();
                int ind = 0;
                for (float l : wordL) {
                    wordinput.push_back(guess.substr(ind, l));
                    ind += l;
                }
                for (string word : wordinput) cout << word << " ";
                cout << endl;
                GetMatchMat(guess, matchmat);
                cout << "compmat " << 1.f / CompareMat(probmat, matchmat) << endl;
                cout << "wordscore " << WordScore(wordinput) << endl;
                cout << "ngramscore " << NGramScore(wordinput) << endl;
                cout << "Acceptance rate is " << (float)acc / prop << endl;
                cout << "Trial stat " << biastrialpos / trials << endl << endl;
            }
        }
    }
}

static map<float, string> solutions;
static float best = -1;
static int totalIters = 0;
void RecursiveInfer(vector<int>& maskKey, vector<string>& answers, vector<string>& words, 
                            vector<vector<string>>& wordGuesses, int level) {
    ++totalIters;
    if (!(totalIters % 1000))
        cout << totalIters << endl;
    // Find smallest list in wordGuesses that isn't 1. If none, return!!
    int easyInd, poss = INT_MAX;
    bool done = true;
    for (int i = 0; i < wordGuesses.size(); ++i) {
        //if (answers[i].empty() && words[i].length() < poss) {
        if (answers[i].empty() && wordGuesses[i].size() < poss) {
            //poss = words[i].length();
            poss = wordGuesses[i].size();
            easyInd = i; 
        }
        if (answers[i].empty()) done = false;
        if (wordGuesses[i].size() == 0) {
            //cout << "backing out" << endl;
            return;
        }
    }
    if (done) {
        float score = WordScore(answers) * 5.f + NGramScore(answers);
        if (score >= 0) {
            string sol = "";
            for (string &s : answers)
                sol += s + " ";
            solutions[score] = sol;
            if (score > best) {
                cout << sol << ": " << score << endl;

                for (string &s : answers) {
                    cout << s << " " << freqs[s] << " ";
                }
                cout << endl;
            }
            if (score > best) best = score;

        }
        if (done) return;
    }
    // Loop o'er all strings in smallest list, recursively calling self
    for (int g = 0; g < wordGuesses[easyInd].size(); ++g) {
//    for (string &guess : wordGuesses[easyInd]) {
        string &guess = wordGuesses[easyInd][g];
        vector<int> maskKeyCopy(maskKey);
        vector<string> answersCopy(answers);
        answersCopy[easyInd] = guess;
        vector<vector<string>> wordGuessesCopy(wordGuesses);
        for (int i = 0; i < guess.size(); ++i)
            maskKeyCopy[words[easyInd][i] - 'a'] = guess[i] - 'a';
        for (int i = 0; i < words.size(); ++i) {
            string mask = PermuteString(maskKeyCopy, words[i]);
            string range = "";
            for (char c : maskKeyCopy)
                if (c != '.' - 'a')
                    range += (char)(c + 'a');
            //wordGuessesCopy[i] = MaskCharsWithBanned(mask, range, wordGuessesCopy[i]);
            wordGuessesCopy[i] = MaskChars(mask, wordGuessesCopy[i]);
        }
        RecursiveInfer(maskKeyCopy, answersCopy, words, wordGuessesCopy, level + 1);
    }
}

void GetCompatibleWords(vector<string> &words, vector<vector<string>> &wordGuesses) {
    for (string &s : words) {
        wordGuesses.push_back(vector<string>());
        for (auto &p : freqs) {
            string word = p.first;
            if (word.length() != s.length() || p.second < FREQ_CUTOFF) continue;
            map<char, char> key;
            bool ok = true;
            for (int i = 0; i < s.length() && ok; ++i) {
                if (key.find(s[i]) == key.end())
                    key[s[i]] = word[i];
                ok = (key[s[i]] == word[i]);
            }
            if (ok) wordGuesses.back().push_back(word);
        }
    }
}

// once upon a time a really amazing thing happened in the capital of this country
void Infer(string input) {

    stringstream ss(input);
    vector<string> words;
    while (ss.good()) {
        string next;
        ss >> next;
        cout << next << ": " << freqs[next] << ", ";
        if (freqs[next] < FREQ_CUTOFF) {
            cout << "WARNING: " << next << " only appears " << freqs[next] << endl << endl;
            exit(0);
        }
        words.push_back(next);
    }
    cout << endl;
    vector<int> key;
    AlphaPermute(key);
    PermuteStringVector(key, words);
    for (string &s : words) cout << s << ' ';
    cout << endl;

    vector<vector<string>> wordGuesses;
    GetCompatibleWords(words, wordGuesses);

    // for (string &s : words)
    //     wordGuesses.push_back(cipherCodes[CipherCode(s)]);

    for (auto &v : wordGuesses) {
        cout << v.size() << ' ';
        sort(v.begin(), v.end(), [] (const string &a, const string &b) {
            return freqs[a] > freqs[b];
        });
    }
    cout << endl;
    vector<int> maskKey(26, '.' - 'a');
    vector<string> answers(words.size(), "");

    RecursiveInfer(maskKey, answers, words, wordGuesses, 0);
    cout << "Total iterations: " << totalIters << endl;
    cout << "SORTED ORDER: " << endl;

    auto iter = solutions.end();
    for (int i = 0; i < 10 && i < solutions.size(); ++i, --iter) ;
    while (iter != solutions.end())
        cout << iter->second << ": " << (iter++)->first << endl;
    /*
    for (auto &s : solutions)
        cout << s.second << ": " << s.first << endl;
        */

}

int main(int argc, const char *argv[]) 
{
    const int n = 3;
    srand (time(NULL));
    //string prefix(argv[1]);
    string prefix = "corpus";
    prefix += "_";
    string freqfile = prefix + "freqs";
    string gramfile = prefix + "3grams";
    //ReadStringIntMap("LiveWrite/COCA_5000.txt", freqs);

    string dictfile = "wordsEn.txt";
    ReadStringIntMap(dictfile, freqs, false);
    map<string, int> corpus;
    ReadStringIntMap(freqfile, corpus);
    ReadStringIntMap(gramfile, ngrams);
    for (auto &p : freqs)
            p.second += corpus[p.first];
    /*
    for (auto &p : freqs) {
        if (dict.find(p.first) == dict.end()) {
            cout << p.first << " ain't in the dictionary!! " << endl;
            ++d;
        } else ++both;
    }
    cout << "SWITCHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    for (auto &p : dict) {
        if (freqs.find(p.first) == freqs.end()) {
            cout << p.first << " ain't in the corpus!! " << endl;
            ++c;
        }
    }
    cout << both << " in both" << endl;
    cout << d << " in corpus and not dict, out of " << freqs.size() << endl;
    cout << c << " in dict and not corpus, out of " << dict.size() << endl;
    return 0;*/

    for (auto &p : freqs) 
        cipherCodes[CipherCode(p.first)].push_back(p.first);

    cout << "Mapped cipher codes" << endl;

    cout << "Input a one line phrase: " << endl;
    string inputline;
    getline(cin, inputline);

    Infer(inputline);
    return 0;


    string inputnos = "";
    vector<float> wordL;
    wordL.push_back(0);
    for (char c : inputline) {
        if (c == ' ') {
            wordL.push_back(0);
        } else {
            inputnos += c;
            wordL.back()++;
        }
    }
    int nchar = inputnos.size();
    vector<float> probmat(nchar * nchar);
    for (int j = 0; j < nchar; ++j) {
        for (int i = 0; i < nchar; ++i) {
            float p = 1.f;
            if (j > i) p = probmat[i * nchar + j];
            if (i > j) {
                float match = (inputnos[i] == inputnos[j] ? 1.f : 0.f);
                float bias = 0.f;
                p = match * (1.f - bias) + randuni() * bias;
                //if (trial(.1f) && p == 0.f) p = randuni() * .5f + .5f;
            }
            probmat[j * nchar + i] = p;
            //cout << inputnos[i] << inputnos[j] << p << "  ";
        }
        //cout << endl;
    }
    cout << "Input: " << inputnos << endl;
    ProbInf(probmat, nchar, wordL);


    // FILE *fp = fopen(argv[1], "r");
    // if (!fp) {
    // 	printf("Cannot open file \"%s\"\n", argv[1]);
    // 	return 0;
    // }
    // map<string, int> freqs;
    // FillInCorpus(fp, freqs);
    // fclose(fp);
    // cout << "Found " << freqs.size() << " words" << endl;

    // int total = 0;
    // for (auto &p : freqs) total += p.second;
    // cout << "Found " << total << " entries" << endl;
    // map<string, int> ngrams;
    // GetNGrams(freqs, ngrams, n);
    // ofstream outfile(gramfile);
    // for (auto &p : ngrams)
    //     outfile << p.first << " " << p.second << endl;
    // outfile.close();

    //RunInference(50, 5000, freqs, ngrams, input, n);
    return 0;
}
