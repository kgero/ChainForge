import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Handle } from 'react-flow-renderer';
import useStore from './store';
import NodeLabel from './NodeLabelComponent'
import {BASE_URL} from './store';
import { Container, Grid, Select, Radio, NumberInput, Group, Paper, Text, Accordion, SegmentedControl } from '@mantine/core';
import './output-grid.css';

import embeddingsFilePath from './nlp_models/glove.6B.50d.10k.txt';

import winkNLP from 'wink-nlp'
import winkModel from 'wink-eng-lite-web-model'
import BM25Vectorizer from 'wink-nlp/utilities/bm25-vectorizer'
import similarity from 'wink-nlp/utilities/similarity.js'

const nlp = winkNLP(winkModel);
const its = nlp.its;
const as = nlp.as;


// from https://github.com/ianarawjo/ChainForge/blob/main/chainforge/react-server/src/store.js
/** The color palette used for displaying info about different LLMs. */
const llmColorPaletteSaturated = ['#44d044', '#f1b933', '#e46161', '#8888f9', '#33bef0', '#bb55f9', '#f7ee45', '#f955cd', '#26e080', '#2654e0', '#7d8191', '#bea5d1'];
// https://coolors.co/44d044-f1b933-e46161-8888f9-33bef0-bb55f9-f7ee45-f955cd
// +60 brightness
// https://coolors.co/b4ecb4-f9e2ad-f4c0c0-cfcffc-ade5f9-e3bafd-fcf8b5-fdbaeb-a8f3cc-a8baf3-cbcdd3-e5dbec
const llmColorPalette = ['#b4ecb4', '#f9e2ad', '#cfcffc', '#ade5f9', '#e3bafd', '#fcf8b5', '#fdbaeb', '#fdbaeb', '#a8f3cc', '#a8baf3', '#cbcdd3', '#e5dbec'];
// +60 brightness again
// https://coolors.co/e2f8e2-fdf4df-fae5e5-ececfe-dff5fd-f4e3fe-fefce2-fee3f7-ddfaeb-dde4fa-eaebed-f4f0f7
const llmColorPaletteDesaturated = ['#e2f8e2', '#fdf4df', '#fae5e5', '#ececfe', '#dff5fd', '#f4e3fe', '#fefce2', '#fee3f7', '#ddfaeb', '#dde4fa', '#eaebed', '#f4f0f7'];

/** The color palette used for displaying variations of prompts and prompt variables (non-LLM differences).
 * Distinct from the LLM color palette in order to avoid confusion around what the data means.
 * Palette adapted from https://lospec.com/palette-list/sness by Space Sandwich */
const varColorPalette = ['#0bdb52', '#e71861', '#7161de', '#f6d714', '#80bedb', '#ffa995', '#a9b399', '#dc6f0f', '#8d022e', '#138e7d', '#c6924f', '#885818', '#616b6d'];

// inspired by https://github.com/sauravjoshi23/GloVe-Embeddings-NLP-JS/tree/master
// returns an object where each key is a word and the value is the vector
const load_embeddings = (textdata) => {
    // let data = fs.readFileSync('./nlp_models/glove.6B.50d.10k.txt')
    let data = textdata;
    data = data.toString()
    data = data.split("\n")

    let glove_embeddings = {}
    for(var i=0; i<data.length;i++){
        let values = data[i]
        //tokenize -> separate word and vector
        let tokens = values.split(' ')
        let word = tokens.slice(0,1)[0];
        let vect = tokens.slice(1,);
        // coverting string to int
        let int_vect = []
        for(var j=0; j<vect.length;j++){
            var num = Number(vect[j]);
            int_vect.push(num);
        }
        glove_embeddings[word] = int_vect;
    }
    return glove_embeddings
}


const dot_product = (vector1, vector2) => {
    if (vector1 == null || vector2 == null) {
        return null;
    }
    let result = 0;
    let length = vector1.length;
    for (let i = 0; i < length; i++) {
      result += vector1[i] * vector2[i];
    }
    return result;
}

const den_calculation = (vector) => {
    let sum = 0;
    let length = vector.length;
    for (let i = 0; i < length; i++) {
      sum += vector[i] * vector[i];
    }
    let result = Math.sqrt(sum);
    return result;
}

const cosine_similarity = (vector1, vector2) => {
    let num = dot_product(vector1, vector2);
    if (num == null) {return null;}
    let den1 = den_calculation(vector1);
    let den2 = den_calculation(vector2);
    let den = den1*den2;
    let cos_sim = num/den;
    return cos_sim;
}


const stop = [",", ";", ":", "?", "!", "(", ")", "i", "me", "my", "myself", "we", "our", "ours", "ourselves", "you", "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself", "she", "her", "hers", "herself", "it", "its", "itself", "they", "them", "their", "theirs", "themselves", "what", "which", "who", "whom", "this", "that", "these", "those", "am", "is", "are", "was", "were", "be", "been", "being", "have", "has", "had", "having", "do", "does", "did", "doing", "a", "an", "the", "and", "but", "if", "or", "because", "as", "until", "while", "of", "at", "by", "for", "with", "about", "against", "between", "into", "through", "during", "before", "after", "above", "below", "to", "from", "up", "down", "in", "out", "on", "off", "over", "under", "again", "further", "then", "once", "here", "there", "when", "where", "why", "how", "all", "any", "both", "each", "few", "more", "most", "other", "some", "such", "no", "nor", "not", "only", "own", "same", "so", "than", "too", "very", "s", "t", "can", "will", "just", "don", "should", "now"];


const getColorFromPalette = (palette, index) => {
    if (index < palette.length) {
        return palette[index];
    }
    return palette[index % palette.length];
}

const GridInspectNode = ({ data, id }) => {

  let is_fetching = false;

  const [embeddings, setEmbeddings] = useState('');

  useEffect(() => {
    fetch(embeddingsFilePath)
      .then(response => response.text())
      .then(data => {
        let embdata = load_embeddings(data);
        setEmbeddings(embdata);
      });
  }, []);  // Empty dependency array

  const [visualization, setVisualization] = useState([]);
  const [jsonResponses, setJSONResponses] = useState(null);

  const [pastInputs, setPastInputs] = useState([]);
  const inputEdgesForNode = useStore((state) => state.inputEdgesForNode);
  const setDataPropsForNode = useStore((state) => state.setDataPropsForNode);

  const [columnValue, setColumnValue] = useState([]);
  const [rowValue, setRowValue] = useState([]);
  const [altValues, setAltValues] = useState([]);
  const [highlightRadioValue, setHighlightRadioValue] = useState([]);
  const [groupColor, setGroupColor] = useState([]);
  const [segmentedViewValue, setSegmentedViewValue] = useState('grid');



  // Update the visualization whenever 'jsonResponses' changes:
  useEffect(() => {
    if (!jsonResponses || (Array.isArray(jsonResponses) && jsonResponses.length === 0))
        return;

    // === Construct a visualization using jsonResponses here ===
    // ....


    console.log("similarity of dog and cat", cosine_similarity(embeddings["cat"], embeddings["dog"]));
    console.log("similarity of dog and table", cosine_similarity(embeddings["table"], embeddings["dog"]));


     // Returns length of longest common substring of X[0..m-1] and Y[0..n-1]
    const LCSubStr = (X, Y) => {
        // Create a table to store lengths of longest common suffixes of substrings.
        // Note that LCSuff[i][j] contains length of longest common suffix of
        // X[0..i-1] and Y[0..j-1]. The first row and first column entries have no
        // logical meaning, they are used only for simplicity of program
        // based on https://www.geeksforgeeks.org/longest-common-substring-dp-29/

        const m = X.length;
        const n = Y.length;

        var LCSuff = Array(m + 1).fill().map(()=>Array(n + 1).fill(0));

        // To store length of the longest common substring
        var result = 0;

        // To store the index of the cell which contains the maximum value.
        // This cell's index helps in building up the longest common substring
        // from right to left.
        let row = 0, col = 0;

        // Following steps build LCSuff[m+1][n+1] in bottom up fashion
        for (let i = 0; i <= m; i++) {
            for (let j = 0; j <= n; j++) {
                if (i == 0 || j == 0)
                    LCSuff[i][j] = 0;
                else if (X[i - 1] == Y[j - 1]) {
                    LCSuff[i][j] = LCSuff[i - 1][j - 1] + 1;
                    if (result < LCSuff[i][j]) {
                        result = LCSuff[i][j];
                        row = i;
                        col = j;
                    }
                } else
                    LCSuff[i][j] = 0;
            }
        }
        if (result == 0) {
            console.log("No common substring."); // DEAL WITH THIS CASE
        }
        let resultStr = [];
        while (LCSuff[row][col] != 0) {
            resultStr.unshift(X[row-1]);
            --result;
            row--;
            col--;
        }
        return {"string": resultStr, "x": row, "y": col};
    }

    /*
       responses is an array of responseTokenized;
       k is the number of unique text strings to return;

       returns dict of overlapping strings and their locations in the format:

       {"text of string": [
            {i: index of response string occurs in,
             j: index of token in response where string starts,
             n: number of tokens in string}
       ]}
    */
    const findOverlap = (responses, k) => {
        if (responses.length <= 1) { return {}; }
        let overlaps = {}
        for (let i = 0; i < responses.length-1; i++) {
            for (let j = i+1; j < responses.length; j++) {
                let lcs = LCSubStr(responses[i], responses[j]);
                // check if lcs crosses sentence boundary
                // either return final object or two final objects to add
                const crossSentBoundary = (respIndex, tokenIndex, numTokens, lcsString) => {
                    // console.log("crossSentBoundary", respIndex, tokenIndex, numTokens, lcsString);
                    const sentObj = jsonResponsesMod[respIndex].sentences;
                    // console.log("sentObj", sentObj);
                    for (let k=0; k < sentObj.length; k++) {
                        // console.log("sentObj[k]", sentObj[k]);
                        // console.log("sentObj[k+1]", sentObj[k+1]);
                        if (!sentObj[k].tokenIndices.includes(tokenIndex)) {
                            // string is not in this sentence
                            continue;
                        } else if (sentObj[k].tokenIndices.includes(tokenIndex+numTokens-1)) {
                            // doesn't cross sentence boundary
                            return [{"i": respIndex, "j": tokenIndex, "n": numTokens, "str": lcsString}];
                        } else {
                            // crosses sentence boundary
                            const firstN = sentObj[k].tokenIndices.length - sentObj[k].tokenIndices.indexOf(tokenIndex);
                            const secondN = numTokens - firstN;
                            const secondTokenIndex = sentObj[k+1].tokenIndices[0];
                            const firstString = jsonResponsesMod[respIndex].responseTokenized.slice(tokenIndex, tokenIndex + firstN);
                            const secondString = jsonResponsesMod[respIndex].responseTokenized.slice(secondTokenIndex, secondTokenIndex + secondN);
                            // console.log('crossSentBoundary', k, lcsString, firstN, secondN);
                            return [
                                {"i": respIndex, "j": tokenIndex, "n": firstN, "str": firstString},
                                {"i": respIndex, "j": secondTokenIndex, "n": secondN, "str": secondString}
                                ];
                        }
                    }

                }
                const firstAdd = crossSentBoundary(i, lcs.x, lcs.string.length, lcs.string);
                const secondAdd = crossSentBoundary(j, lcs.y, lcs.string.length, lcs.string);
                const adds = firstAdd.concat(secondAdd);
                // console.log("adds", adds);

                for (let a=0; a<adds.length; a++) {
                    if (!overlaps.hasOwnProperty(adds[a].str)) {
                        overlaps[adds[a].str] = [];
                    }
                    overlaps[adds[a].str].push(adds[a])
                }

                // Old implementation before breaking on sentence boundaries

                // if (!overlaps.hasOwnProperty(lcs.string)) {
                //     overlaps[lcs.string] = []
                // }
                // overlaps[lcs.string].push({"i": i, "j": lcs.x, "n": lcs.string.length})
                // overlaps[lcs.string].push({"i": j, "j": lcs.y, "n": lcs.string.length})
                // overlaps[lcs.string]
            }
        }
        for (const key in overlaps) {
            if (overlaps[key].length === 1) {
                delete overlaps[key];
            }
        }
        console.log("overlaps", overlaps);
        // Convert to [key, value] pairs and calculate max n for each overlap
        let pairs = Object.entries(overlaps).map(([string, locations]) => {
            let maxN = Math.max(...locations.map((location) => location.n));
            let count = locations.length;
            return [string, locations, maxN, count];
        }).filter(([, , maxN]) => maxN > 2);
        // console.log("pairs", pairs);

        // Remove overlaps that have n < 3

        // Sort by max n, in descending order, then by count, also in descending order
        // !! reversing this, sort first by count, then by max n
        // pairs.sort((a, b) => b[2] - a[2] || b[3] - a[3]);
        // pairs.sort((a, b) => b[3] - a[3] || b[2] - a[2]);

        // try a weighted combo, instead of ordering by one and then the other
        // suggested by ChatGPT lol
        const weightForMaxN = 0.5;  // You can adjust this value to emphasize or de-emphasize maxN
        const weightForCount = 1 - weightForMaxN;  // Makes sure the two weights sum up to 1

        pairs.sort((a, b) =>
            (b[2] * weightForMaxN + b[3] * weightForCount) -
            (a[2] * weightForMaxN + a[3] * weightForCount)
        );

        // Get top k overlaps
        let topKoverlaps = {};
        let blockedIndices = Array(responses.length)
            .fill()
            .map(() => new Set());

        for (let i = 0; i < pairs.length; ++i) {
            const [string, locations] = pairs[i];

            if (locations.length === 1) {
              continue;
            }

            // check if a string is overlapping with an already-added string
            let overlapped = false;

            for (const location of locations) {
              for (let idx = location.j; idx < location.j + location.n; idx++) {
                if (blockedIndices[location.i].has(idx)) {
                  overlapped = true;
                  break;
                }
              }
              if (overlapped) break;
            }

            if (!overlapped) {
              topKoverlaps[string] = locations;
              if (Object.keys(topKoverlaps).length === k) break;

              for (const location of locations) {
                for (let idx = location.j; idx < location.j + location.n; idx++) {
                  blockedIndices[location.i].add(idx);
                }
              }
            }
        }


        return topKoverlaps;
    }

    // return a new json object where the keys are the same as the input
    // json and the values are unique hex colors based on llmColorPalette
    const getColors = (json) => {
        let colors = {};
        let counter = 0;
        for (const key in json) {
            colors[key] = llmColorPalette[counter];
            counter++;
        }
        return colors;
    }

    // white space based tokenization
    const tokenize = (string) => {
        string = string.replaceAll("\n", " \n ");
        return string.split(" ");
    }

    console.log("CUSTOM INSPECT NODE");
    console.log('jsonResponses', jsonResponses);
    // first create an object for each response, copying all other attributes
    // (originally 'response' is a list of responses, but creating one per object makes some other logic simpler)
    // also add 'llm' attribute to 'vars' (because we're going to use 'vars' to set the cols and rows)
    // and create a tokenized version of the response
    let counterResp = 0;
    let jsonResponsesMod = jsonResponses.flatMap(obj =>
        obj.responses.map((response, index) => {
            let cleanresponse = response.replaceAll("\n\n", " \n\n<br/><br/> ");
            // const cleanresponse = response;
            const doc = nlp.readDoc(cleanresponse);
            const sentencesOut = doc.sentences().out();
            const sentenceTokens = sentencesOut.map((text) => tokenize(text));
            const flatTokens = sentenceTokens.flat();
            const lemmas = flatTokens.map((token) => nlp.readDoc(token).tokens().out(its.lemma)[0]);
            let indices = [];
            let counter = 0;

            for (let sentence of sentenceTokens) {
                let sentenceIndices = [];
                for (let token of sentence) {
                    sentenceIndices.push(counter);
                    counter++;
                }
                indices.push(sentenceIndices);
            }

            const sentences = [...Array(sentencesOut.length).keys()].map((i) => ({
                text: sentencesOut[i],
                tokenIndices: indices[i],
                bow: nlp.readDoc(sentencesOut[i].replaceAll("<br/>", "")).tokens().out(its.value, as.bow)
            }));
            let currObj = {
                ...obj,
                vars: { ...obj.vars, responseNum: index.toString(), model: obj.llm },
                response: response,
                responseTokenized: flatTokens,
                responseLemmas: lemmas,
                index: counterResp,
                sentences: sentences
            };
            counterResp++;
            return currObj;
        })
    );
    console.log('jsonResponsesMod', jsonResponsesMod);

    // function to learn the bm25 vectors
    const getBm25Vectors = () => {
        const bm25 = BM25Vectorizer();
        const corpus = jsonResponsesMod.map((obj) => obj.response);
        corpus.forEach((doc) => bm25.learn(nlp.readDoc(doc).tokens().out(its.lemma)));;
        return bm25;
    }

    // calcualte idf arrays for each vector
    const getIdfArray = (bm25) => {
        // pull out tf-idf and top terms
        const idf_array = bm25.out(its.idf);
        const idf_obj = idf_array.reduce((obj, [key, val]) => {
            obj[key] = val;
            return obj;
        }, {});
        return idf_obj;
    }

    // calculate tf-idf for doc i and return top n scoring terms
    const topTfidfTerms = (i, n) => {
        const tf = bm25.doc(i).out(its.tf);
        const tfidf = tf.map((arr) => [arr[0], arr[1]*idf_obj[arr[0]]])
        tfidf.sort((a, b) => b[1] - a[1]);
        let top = tfidf.slice(0, n);
        return top.map(term => term[0]);
    };

    // function to remove stopwords from tfidfArray
    const removeStopWords = (outputArray, stopWords, n) => {
        const updatedArray = outputArray.map(innerArray => {
            return innerArray.filter(word => !stopWords.includes(word));
        });
        return updatedArray;
    }

    const bm25 = getBm25Vectors();
    const idf_obj = getIdfArray(bm25);
    const baseTfidfArray = [...Array(jsonResponsesMod.length).keys()].map((i) => topTfidfTerms(i, 10));
    const tfidfN = 5;
    const tfidfArray = removeStopWords(baseTfidfArray, stop).map((arr) => arr.slice(0,tfidfN));

    console.log('idf_obj', idf_obj);
    console.log('tfidfArray', tfidfArray);


    /*  BETTER CLUSTERING
    */
    const absolute_sim_dist_threshold = 0.55;
    const group_making_threshold = 1.2;
    const a = 1.5; //relative weight on content similarity
    const b = 1; //relative weight on location coherence (less important than content, I think)
    const percentage_same_response_threshold = 0.1;

    // modifying allSentences; doesn't use response_lengths
    const compute_response_lengths_and_add_normalized_locations = () => { 
        let allSentences = jsonResponsesMod.map((respObj) => {
            return respObj.sentences.map((sentObj, idx) => {
                const sentenceDetails = {
                    bow: sentObj.bow,
                    respIndex: respObj.index,
                    sentIndex: idx,
                    sentence: sentObj.text
                }
                return sentenceDetails
            })
        }).flat();

      let response_lengths = {};
      for (let i=0; i<allSentences.length; i++){
        if (Object.keys(response_lengths).includes(allSentences[i].respIndex)) {
          if (response_lengths[allSentences[i].respIndex]<allSentences[i].sentIndex){
            response_lengths[allSentences[i].respIndex] = allSentences[i].sentIndex;
          }
        } else {
          response_lengths[allSentences[i].respIndex]=allSentences[i].sentIndex;
        }
        //return parseFloat(response_lengths[allSentences[i].respIndex])
        
        
      }
      for (let i=0; i<allSentences.length; i++){
        //allSentences[i].normalized_location = parseFloat(allSentences[i].sentIndex)/parseFloat(response_lengths[allSentences[i].respIndex]);
          allSentences[i].total_resp_length = response_lengths[allSentences[i].respIndex];
          if (allSentences[i].sentIndex == 0) {allSentences[i].normalized_location = 0} else {
            allSentences[i].normalized_location = parseFloat(allSentences[i].sentIndex)/parseFloat(allSentences[i].total_resp_length);
          }
          
      }
      return allSentences;
    }

    const distance = (sentence_obj1, sentence_obj2) => {
      var mismatch_total = 0;
      var total_tokens = 0;
      var mismatch_dict = {};
      var match_dict = {};
      // could replace w a different distance measure, could be roughly BoW cosine
      // this could be based on greying
      for (const [key, value] of Object.entries(sentence_obj1.bow)) {
        let value2 = sentence_obj2.bow[key];
        if (value2 === undefined) {
          mismatch_total += value;
          mismatch_dict[key] = value;
          total_tokens += value;
        } else {
          match_dict[key] = {value, value2};
          total_tokens += value + value2;
        }
        
      }
      for (const [key, value] of Object.entries(sentence_obj2.bow)) {
        if (sentence_obj1.bow[key] === undefined) {
          mismatch_total += value;
          mismatch_dict[key] = value;
          total_tokens += value;
        }
      }
      var mismatch_score = mismatch_total/total_tokens;
      return { mismatch_score, mismatch_total, total_tokens, mismatch_dict, match_dict };
    }

    const get_pairs = (minSimScoreThreshold, allSentences) => {

      // then make a similarity matrix
      let simMatrix = Array(allSentences.length)
        .fill(0)
        .map(() => Array(allSentences.length).fill(0)); // TODO: SAFEST THING IS TO FILL WITH INF NOT ZEROS
      let simMatrix_unmodified = Array(allSentences.length)
        .fill(0)
        .map(() => Array(allSentences.length).fill(0)); // TODO: SAFEST THING IS TO FILL WITH INF NOT ZEROS
      for (let i = 0; i < allSentences.length; i++) {
        for (let j = 0; j < i; j++) { //instead of j<i, I'm filling out the whole matrix for the unmodified version (not used for clustering) so it's bidirectional
          simMatrix[i][j] = distance(allSentences[i], allSentences[j])[
            "mismatch_score"
          ];
        }
      }
      for (let i = 0; i < allSentences.length; i++) {
        for (let j = 0; j < allSentences.length; j++) { //instead of j<i, I'm filling out the whole matrix for the unmodified version (not used for clustering) so it's bidirectional
          simMatrix_unmodified[i][j] = distance(allSentences[i], allSentences[j])[
            "mismatch_score"
          ];
        }
      }
      // console.log("simMatrix", simMatrix);

      // run clustering algo
      let clusters = [];
      let cluster_dists = []
      const getMaxSimIndices = (arr) => {
        let min = Infinity;
        let indices = [];
        for (let i = 0; i < arr.length; i++) {
          for (let j = 0; j < i; j++) {   //TODO: EXPAND TO ARR.LENGTH IF ASSYMETRIC FUNCTION
            if (arr[i][j] < min) {
              min = arr[i][j];
              indices = [i, j];
            }
          }
        }
        return [min, indices];
      };
      let minSimScore; let bestPair;
      [minSimScore, bestPair] = getMaxSimIndices(simMatrix);
      clusters.push(bestPair);
      cluster_dists.push(minSimScore);
      simMatrix[bestPair[0]][bestPair[1]] = Infinity;
      
      while (minSimScore<minSimScoreThreshold){
        [minSimScore, bestPair] = getMaxSimIndices(simMatrix);
        clusters.push(bestPair);
        cluster_dists.push(minSimScore);
        simMatrix[bestPair[0]][bestPair[1]] = Infinity;
      }
      
      // only output that is used is simMatrix_unmodified
      return {minSimScore, clusters, cluster_dists, simMatrix_unmodified}
    }

    const getCurrentGroupIdx = (sentenceID,groups) => {
      for (let i = 0; i<groups.length; i++){
        if (groups[i].includes(sentenceID)) {
          return i
        }
      }
      return -1
    }

    const getPercentageSharedResponses = (groupA,groupB) => {
      let responses_included = [];
      for (let i = 0; i<groupA.length; i++){
        responses_included.push(allSentences[groupA[i]].respIndex);
      }
      for (let i = 0; i<groupB.length; i++){
        responses_included.push(allSentences[groupB[i]].respIndex);
      }
      let responses_included_unique = Array.from(new Set(responses_included));
      return 1-(responses_included_unique.length/responses_included.length)
    }

    const addToGroups = (groups,pair) => {
        let group_idx_pair0 = getCurrentGroupIdx(pair[0], groups)
        let group_idx_pair1 = getCurrentGroupIdx(pair[1], groups)
        if (group_idx_pair0 == -1 && group_idx_pair1 == -1) {
          groups.push(pair);
        } else if (group_idx_pair0 == group_idx_pair1){
          return groups //do nothing!
        } else if (group_idx_pair0 > -1 && group_idx_pair1 > -1) {

          //TODO: WHAT PERCENTAGE OF SENTENCES IN THIS CLUSTER SHARE A RESPONSE? IF 30% OR MORE SHARE THE SAME RESPONSE WHEN COMBINED, DON'T COMBINE
          if (getPercentageSharedResponses(groups[group_idx_pair0],groups[group_idx_pair1]) > percentage_same_response_threshold) {
            return groups
          }
          
          //WHICH IS BETTER FOR MINIMIZING DISTANCES OF ADJACENT GROUP MEMBERS---ADDING A TO B OR B TO A?
          let end_of_group_a = groups[group_idx_pair0].slice(-1);
          let start_of_group_b = groups[group_idx_pair1].slice(0,1);
          let end_of_group_b = groups[group_idx_pair1].slice(-1);
          let start_of_group_a = groups[group_idx_pair0].slice(0,1);

          let dist_a_push_b = simMatrix[end_of_group_a][start_of_group_b];
          let dist_b_push_a = simMatrix[end_of_group_b][start_of_group_a];

          if (dist_a_push_b <= dist_b_push_a) {
            let group_part_b = groups[group_idx_pair1];
            for (let i = 0; i<group_part_b.length; i++){
              groups[group_idx_pair0].push(group_part_b[i]);
            }
            groups.splice(group_idx_pair1,1); //removes the group was added elsewhere
          } else {
            let group_part_a = groups[group_idx_pair0];
            for (let i = 0; i<group_part_a.length; i++){
              groups[group_idx_pair1].push(group_part_a[i]);
            }
            groups.splice(group_idx_pair0,1); //removes the group was added elsewhere
          }
          //TODO: CONSIDER TESTING FOR MERGED CLUSTER PROPERTIES/STATS
        } else if (group_idx_pair0 > -1){

          if (getPercentageSharedResponses(groups[group_idx_pair0],[pair[1]]) > percentage_same_response_threshold) {
            return groups
          }
          
          let start_of_group = groups[group_idx_pair0].slice(0,1);
          let end_of_group = groups[group_idx_pair0].slice(-1);
          let dist_at_end_of_group = simMatrix[end_of_group][pair[1]];
          let dist_at_start_of_group = simMatrix[pair[1]][start_of_group];
          if (dist_at_start_of_group <= dist_at_end_of_group){
            groups[group_idx_pair0].splice(0,0,pair[1]);
          } else{
            groups[group_idx_pair0].splice(groups[group_idx_pair0].length,0,pair[1]);
          }
          //TODO: CONSIDER TESTING FOR MERGED CLUSTER PROPERTIES/STATS
          //TODO: WHAT PERCENTAGE OF SENTENCES IN THIS CLUSTER SHARE A RESPONSE? IF 30% OR MORE SHARE THE SAME RESPONSE WHEN COMBINED, DON'T COMBINE
        } else if (group_idx_pair1 > -1){

          if (getPercentageSharedResponses(groups[group_idx_pair1],[pair[0]]) > percentage_same_response_threshold) {
            return groups
          }
          
          let start_of_group = groups[group_idx_pair1].slice(0,1);
          let end_of_group = groups[group_idx_pair1].slice(-1);
          let dist_at_end_of_group = simMatrix[end_of_group][pair[0]];
          let dist_at_start_of_group = simMatrix[pair[0]][start_of_group];
          if (dist_at_start_of_group <= dist_at_end_of_group){
            groups[group_idx_pair1].splice(0,0,pair[0]);
          } else {
            groups[group_idx_pair1].splice(groups[group_idx_pair1].length,0,pair[0]);
          }
          //TODO: CONSIDER TESTING FOR MERGED CLUSTER PROPERTIES/STATS
          //TODO: WHAT PERCENTAGE OF SENTENCES IN THIS CLUSTER SHARE A RESPONSE? IF 30% OR MORE SHARE THE SAME RESPONSE WHEN COMBINED, DON'T COMBINE
        }
      return groups; //TODO: PUT IN LOGIC
    }

    const get_mean = (arr) => {
        const sum = arr.reduce((a, b) => a + b, 0);
        return (sum / arr.length);
    }

    const get_ordered_groups = (absolute_sim_dist_threshold, allSentences) => {
      // this outputs the ordered clusters
      let groups = [];
      let ordered_pairs = get_pairs(absolute_sim_dist_threshold, allSentences).clusters;
      for (let pair_idx = 0; pair_idx < ordered_pairs.length; pair_idx++) {
        let pair = ordered_pairs[pair_idx];
        let pair_loc_abs_diff = Math.abs(
          allSentences[pair[0]].normalized_location -
            allSentences[pair[1]].normalized_location
        );
        let pair_sim_dist = simMatrix[pair[0]][pair[1]]; // TODO: SAFEST IS TO CONSIDER BOTH DIRECTIONS OR MIN OF BOTH DIRECTIONS IF DISTANCE IS ASSYMMETRIC
        if (a * pair_sim_dist + b * pair_loc_abs_diff < group_making_threshold) {
          groups = addToGroups(groups, pair);
          //return groups
        }
      }
      // add singleton clusters to groups
      let flattened_groups = groups.flat();
      for (let i=0; i<allSentences.length; i++){
        if (!flattened_groups.includes(i)){
          groups.push([i]); // adds a singleton
        }
      }
      //return groups

      //now order it by average or median normalized location
      //compute aggregate normalized location for each group
      let ordered_groups = [];
      let median_groups_locations = [];
      let num_groups = groups.length;
      for (let group_idx = 0; group_idx< num_groups; group_idx++){
        let group_locations = [];
        for (let intra_group_idx = 0; intra_group_idx<groups[group_idx].length; intra_group_idx++){
          group_locations.push(allSentences[groups[group_idx][intra_group_idx]].normalized_location);
        }
        median_groups_locations.push(get_mean(group_locations));
      }
      //return median_groups_locations
      
      for (let i=0; i<num_groups; i++){
        let next_group_idx = median_groups_locations.indexOf(Math.min(...median_groups_locations)); //TODO: TIE BREAK TO PUT SENTENCES FROM LONGER DOCS FIRST //update: using mean fixed this
        ordered_groups.push(groups[next_group_idx]);
        groups.splice(next_group_idx,1);
        median_groups_locations.splice(next_group_idx,1);
        //return median_groups_locations
      }
      
      return ordered_groups;
    }

    const allSentences = compute_response_lengths_and_add_normalized_locations();
    const simMatrix = get_pairs(absolute_sim_dist_threshold,allSentences).simMatrix_unmodified;
    const pairs = get_pairs(absolute_sim_dist_threshold,allSentences).clusters;
    const clusters = get_ordered_groups(absolute_sim_dist_threshold,allSentences);
    console.log("allSentences", allSentences);
    console.log("clusters", clusters);

    
    // get the original prompt
    let prompt = jsonResponses[0].prompt
    for (const key in jsonResponses[0].vars) {
        prompt = prompt.replace(jsonResponses[0].vars[key], "{"+key+"}");
    }

    // find all vars for table
    const var_names = [...new Set(jsonResponsesMod.map((obj) => Object.keys(obj.vars)).flat(1))];
    const tab_options = [var_names].flat(1);
    console.log("var_names", var_names, "tab_options", tab_options);

    // set up the alt options -- the values that will not be represented in the table axes
    // right now this is same as tab options, but later we'll remove the selected table axes as options
    let alt_options = tab_options.slice();

    // set the var for rows and columns
    // note that columnValue and rowValue are set by a dropdown (see e.g. handleColumnValueChange)
    let col_names = [...new Set(jsonResponsesMod.map((obj) => obj.vars[columnValue]))];
    let row_names = [...new Set(jsonResponsesMod.map((obj) => obj.vars[rowValue]))];

    // look for the data for each cell
    const gridResponses = [];
    for (let y=0; y<row_names.length; y++) {
        let row = [];
        for (let x=0; x<col_names.length; x++) {
            // find the data that goes in that cell given the row and col
            // assumes it should only find one (i.e. takes first item in filtered list)
            let cellData = jsonResponsesMod.filter((resp) =>
                resp.vars[columnValue] === col_names[x]
                && resp.vars[rowValue] === row_names[y]);
            // now filter it given any extra options that need to be selected
            cellData = cellData.filter(resp => {
                return Object.entries(altValues).every(([key, val]) => {
                    if (key === columnValue) { return resp; }
                    if (key === rowValue) { return resp; }
                    return resp.vars[key] === val;
                });
            });
            cellData = cellData[0];
            if (!cellData) {
                cellData = {'responseTokenized': ['∅']};
            }
            row.push(cellData);
        }
        gridResponses.push(row);
    }
    console.log('gridResponses', gridResponses);
    console.log('columnValue', columnValue);
    console.log('altValues', altValues);

    // calculate longest common substring for each pair
    // const rowLCS = gridResponses.map((row) => findOverlap(row.map((cell) => cell.responseTokenized), 5));
    // const colLCS = gridResponses[0].map((_, i) => findOverlap(gridResponses.map((row) => row[i].responseTokenized), 5));
    // const allLCS = findOverlap(gridResponses.flatMap(item => item).map(item => item.responseTokenized), 5);
    const numOverlaps = Math.min(Math.floor(jsonResponsesMod.length / 2), llmColorPalette.length-1);
    const allLCS = findOverlap(jsonResponsesMod.map(item => item.responseTokenized), numOverlaps);
    const allLCS_colors = getColors(allLCS);
    console.log("allLCS", allLCS);
    console.log("allLCS_colors", allLCS_colors);


    // lcs highlighting uses the same functions as the old ngram ones
    // because the detecting algorithms return the same output form
    // const shouldHighlightRowNgrams = (rowIndex, cellIndex, tokenIndex) => {
    //     for (const [ngram, locations] of Object.entries(rowLCS[rowIndex])) {
    //         for (let l = 0; l < locations.length; l++) {
    //           if (locations[l].i === cellIndex) {
    //             // if this ngram occurs in this string
    //             if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
    //                 return true;
    //             }
    //           }
    //         }
    //     }
    //     return false;
    // }
    // const shouldHighlightColNgrams = (rowIndex, cellIndex, tokenIndex) => {
    //     for (const [ngram, locations] of Object.entries(colLCS[cellIndex])) {
    //         for (let l = 0; l < locations.length; l++) {
    //           if (locations[l].i === rowIndex) {
    //             // if this ngram occurs in this string
    //             if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
    //                 return true;
    //             }
    //           }
    //         }
    //     }
    //     return false;
    // }
    const shouldHighlightAllNgrams = (cell, tokenIndex) => {
        for (const [ngram, locations] of Object.entries(allLCS)) {
            for (let l=0; l<locations.length; l++) {
                if (locations[l].i === cell.index) {
                    if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
                        return [true, allLCS_colors[ngram]];
                    }
                }
            }
        }
        return [false, null];
    }
    const shouldHighlightTfidf = (cell, tokenIndex) => {
        const currTokenLemma = cell.responseLemmas[tokenIndex];
        if (tfidfArray[cell.index].includes(currTokenLemma)) {
            let colorIndex = cell.index;
            // if (colorIndex >= llmColorPalette.length) {
            //     colorIndex = colorIndex - llmColorPalette.length;
            // }
            // return [true, llmColorPalette[colorIndex]];
            return [true, "#FFFF8F"];
        }
        return [false, null];
    }
    const shouldHighlightSentence = (cell, tokenIndex) => {
        if (cell.sentences == undefined) {
            return [false, null];
        }

        // first get the sentence index of the token
        let sentIndex;
        for (let i=0; i<cell.sentences.length; i++) {
            if (cell.sentences[i].tokenIndices.includes(tokenIndex)) {
                sentIndex = i;
                break;
            }
        }
        // then get the sentence id from allSentences
        let sentId;
        for (let i=0; i<allSentences.length; i++) {
            if (allSentences[i].respIndex == cell.index && allSentences[i].sentIndex == sentIndex) {
                sentId = i;
                break;
            }
        }

        // then check if sentId is in a cluster
        let clusterId = null;
        for (let i=0; i<clusters.length; i++) {
            if (clusters[i].includes(sentId)) {
                clusterId = i;
                break;
            }
        }


        if (clusterId !== null && clusters[clusterId].length > 1) {
            return [true, getColorFromPalette(llmColorPalette, clusterId)]
        }

        return [false, null];
    }


    col_names.unshift(""); // add empty string for row headers

    const header = col_names.map((el) => { return <th>{el}</th> });

    const cells = gridResponses.map((row, rowIndex) => {
        return (
            <tr>
                <th>{row_names[rowIndex]}</th>
                {row.map((cell, cellIndex) => (
                    <td key={cellIndex}>
                            {cell.responseTokenized.map((token, tokenIndex) => {
                                let highlight = false; 
                                let oddEven = false;
                                let color = null;
                                let spanStyle;

                                if (highlightRadioValue === "lcs") {
                                    [highlight, color] = shouldHighlightAllNgrams(cell, tokenIndex);
                                }
                                else if (highlightRadioValue === "tfidf") {
                                    [highlight, color] = shouldHighlightTfidf(cell, tokenIndex);
                                    // color = null;
                                    // if (highlight) {
                                    //     spanStyle = {fontWeight: "bold"}
                                    //     highlight = false;
                                    //     console.log("tfidf highlighting");
                                    // }
                                    
                                }
                                else if (highlightRadioValue === "sent") {
                                    [highlight, color] = shouldHighlightSentence(cell, tokenIndex);
                                }

                                // spanStyle = highlight ? {backgroundColor: oddEven ? 'thistle' : 'plum'} : {};
                                
                                if (color) {
                                    spanStyle = {backgroundColor: color}
                                }

                                if (token == "<br/>") {return <br/>;}
                                if (token == "<br/><br/>") {return <><br/><br/></>;}
                                return <span key={tokenIndex} style={spanStyle}>{token} </span>;;
                            })}
                    </td>
                ))}
            </tr>
        );
    });

    const table = (<table class='outputgrid'>
            <tbody>
                <tr>{header}</tr>
                {cells}
            </tbody>
        </table>);

    const getSentObj = (sentId) => {
        const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
        const sentIndex = allSentences[sentId].sentIndex;
        return respObj.sentences[sentIndex];
    }
    const getToken = (sentId, tokenIndex) => {
        const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
        return respObj.responseTokenized[tokenIndex];
    }

    const groupColorOptions = [...new Set(jsonResponsesMod.map((obj) => obj.vars[groupColor]))];

    const groupColorOptionsDisplay = groupColorOptions.map((colorOption, colorIndex) => {
        const spanStyle = {backgroundColor: getColorFromPalette(llmColorPalette, colorIndex)};
        const thisKey = {colorOption} + "-" + {colorIndex}
        return (
            <>
                <span style={spanStyle}>&nbsp;&nbsp;&nbsp;</span>
                <span>&nbsp;</span>
                <span key={thisKey}>{colorOption}</span>
                <span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>
            </>)
    })

    const countGreyedWords = (sentId1, sentId2) => {

        const getSentObj = (sentId) => {
            const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
            const sentIndex = allSentences[sentId].sentIndex;
            return respObj.sentences[sentIndex];
        }
        const getToken = (sentId, tokenIndex) => {
            const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
            return respObj.responseTokenized[tokenIndex];
        }

        // set up all the relevant object properties
        const respObj = jsonResponsesMod[allSentences[sentId1].respIndex];
        const sentIndex = allSentences[sentId1].sentIndex;
        const tokenIndices = respObj.sentences[sentIndex].tokenIndices; // this is the thing we really care about
        const filteredTokenIndices = tokenIndices.filter((tokenIndex) => getToken(sentId1, tokenIndex) !== "<br/><br/>");

        let count = 0;
        
        filteredTokenIndices.map((tokenIndex, mappedTokenIndex) => {
            const token = getToken(sentId1, tokenIndex);

            const prevSentObj = getSentObj(sentId2);
            const prevTokenIndices = prevSentObj.tokenIndices;
            const filteredPrevTokenIndices = prevTokenIndices.filter((tokenIndex) => getToken(sentId2, tokenIndex) !== "<br/><br/>");

            let prevTokenIndex = filteredPrevTokenIndices[mappedTokenIndex];
            let prevToken = getToken(sentId2, prevTokenIndex);

            if (prevToken === token) {
                count++;
            }
        });

        // return count; // for absolute value
        return count/filteredTokenIndices.length; // if you'd prefer a percentage
        // note that you need to decide if you want to normalize on the length
        // of the first or second sentence...
    }

    const grouping = clusters.map((cluster, clusterIndex) => {
        const clusterSentences = cluster.map((sentId, mappedSentIndex, arr) => {
            const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
            const respColorIndex = groupColorOptions.indexOf(respObj.vars[groupColor]);
            const sentIndex = allSentences[sentId].sentIndex;
            const sentObj = respObj.sentences[sentIndex];
            const tokenIndices = respObj.sentences[sentIndex].tokenIndices;
            const filteredTokenIndices = tokenIndices.filter((tokenIndex) => getToken(sentId, tokenIndex) !== "<br/><br/>");

            // for testing function that counts greyed out words
            // if (mappedSentIndex > 0) {
            //     const numGreyedTokens = countGreyedWords(sentId, arr[mappedSentIndex-1]);
            //     console.log('numGreyedTokens', sentId, numGreyedTokens)
            // }
            
            return (
                <p key={sentId} className="sentenceP">
                <span style={{backgroundColor: getColorFromPalette(llmColorPalette, respColorIndex)}}>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>
                <span>&nbsp;</span>
                {filteredTokenIndices.map((tokenIndex, mappedTokenIndex) => {

                    const token = getToken(sentId, tokenIndex);
                    let spanStyle = {}

                    if (mappedSentIndex > 0) {
                        const prevSentId = arr[mappedSentIndex-1];
                        const prevSentObj = getSentObj(prevSentId);
                        const prevTokenIndices = prevSentObj.tokenIndices;
                        const filteredPrevTokenIndices = prevTokenIndices.filter((tokenIndex) => getToken(prevSentId, tokenIndex) !== "<br/><br/>");

                        let prevTokenIndex = filteredPrevTokenIndices[mappedTokenIndex];
                        let prevToken = getToken(prevSentId, prevTokenIndex);

                        if (prevToken === token) {
                            spanStyle = {color: 'grey'}
                        }
                    }
                    if (token == "<br/>") {return <></>;}
                    if (token == "<br/><br/>") {return <></>;}
                    return (<span style={spanStyle}>{token} </span>);
                })}
                </p>
                );
        });
        return <div className="clusterGroup" key={clusterIndex}>{clusterSentences}</div>;
    });


    

    // ==========================================================

    const handleColumnValueChange = (new_val) => {
        setColumnValue(new_val);
        setAltValues({});
    };
    const handleRowValueChange = (new_val) => {
        setRowValue(new_val);
        setAltValues({});
    };
    const handleAltValuesChange = (new_val, name) => {
        setAltValues(prevState => {
            return {...prevState, [name]: new_val}
        });
    };
    const resetAltValues = () => {
        setAltValues({});
      };

    const handleHighlightRadioValue = new_val => {
        setHighlightRadioValue(new_val);
    };

    const handleSegmentedViewValueChange = (new_val) => {
        setSegmentedViewValue(new_val);
    };

    const handleGroupColorChange = (new_val) => {
        setGroupColor(new_val);
    };

    const gridHighlightingOptions = (<div>
        <div style={{fontSize: '10pt', color: '#777'}}>Select what to highlight:</div>
        <Grid>
            <Grid.Col span={10}>
                <Radio.Group
                  value={highlightRadioValue}
                  onChange={handleHighlightRadioValue}
                  name="highlightRadioValue"
                  defaultValue="none"
                >
                  <Group mt="xs">
                    <Radio value="none" label="None" />
                    <Radio value="lcs" label="Exact Matches" />
                    <Radio value="tfidf" label="Unique Words" />
                    <Radio value="sent" label="Similar Sentences" />
                  </Group>
                </Radio.Group>
            </Grid.Col>
        </Grid>
    </div>);

    const groupHighlightingOptions = (<>
        <Grid.Col span={2} style={{backgroundColor: '#eee', borderRadius: '5px', margin: '10px'}}>
        <div style={{fontSize: '10pt', color: '#777'}}>Color responses by:</div>
            <Select clearable
              size="xs"
              onChange={handleGroupColorChange}
              label=" "
              placeholder="Pick one"
              defaultValue="Model"
              data={tab_options}
              value={groupColor}
            />
        </Grid.Col>
        <Grid.Col span={6} style={{backgroundColor: '#eee', borderRadius: '5px', margin: '10px'}}>
            <div style={{fontSize: '10pt', color: '#777'}}>Legend:</div>
            <p></p>
            {groupColorOptionsDisplay}
        </Grid.Col>
    </>);

    const alt_select_obj = tab_options.map((alt_value, index) => {
        if (columnValue == null || rowValue == null || columnValue == "" || rowValue == "") {
            return (<></>);
        }
        if (alt_value == columnValue || alt_value == rowValue) {
            return (<></>);
        }
        // get the actual options for this value
        let these_options = [...new Set(jsonResponsesMod.map((obj) => obj.vars[alt_value]))];
        return (
            <Grid.Col span={3}>
                <Select
                  size="xs"
                  onChange={(value) => handleAltValuesChange(value, alt_value)}
                  label={alt_value}
                  placeholder="Pick one"
                  defaultValue=""
                  data={these_options}
                  value={altValues[alt_value]}
                />
            </Grid.Col>
        );
    });

    const gridDimensionOptions = (<>
        <Grid.Col span={4} style={{backgroundColor: '#eee', borderRadius: '5px', margin: '10px'}}>
            <div style={{fontSize: '10pt', color: '#777'}}>Set table dimensions:</div>
            <Grid>
                <Grid.Col span={6}>
                    <Select clearable
                      size="xs"
                      onChange={handleColumnValueChange}
                      label="Columns"
                      placeholder="Pick one"
                      defaultValue="Model"
                      data={tab_options}
                      value={columnValue}
                    />
                </Grid.Col>
                <Grid.Col span={6}>
                    <Select clearable
                      size="xs"
                      onChange={handleRowValueChange}
                      label="Rows"
                      placeholder="Pick one"
                      defaultValue="Model"
                      data={tab_options}
                      value={rowValue}
                    />
                </Grid.Col>
            </Grid>
        </Grid.Col>
        <Grid.Col span={5} style={{backgroundColor: '#eee', borderRadius: '5px', margin: '10px'}}>
            <div style={{fontSize: '10pt', color: '#777'}}>Set extra dimensions:</div>
            <Grid>
                {alt_select_obj}
            </Grid>
        </Grid.Col>
    </>);




    // Set the HTML / React element
    const my_vis_component = (<Container fluid="true">
        <div style={{position: "sticky", top: "0", backgroundColor: "white"}}>
            <Grid> 
                <Grid.Col span={2} style={{backgroundColor: '#E7F5FF', borderRadius: '5px', margin: '10px'}}>
                    <div style={{fontSize: '10pt', color: '#777'}}>Set display type:</div>
                    <Grid align="flex-end">
                        <Grid.Col>
                        <label class="mantine-InputWrapper-label mantine-Select-label mantine-jkwmhw" for="mantine-iyl5n76na" id="mantine-iyl5n76na-label" style={{opacity: '0'}}>model</label>
                        <SegmentedControl color="blue"
                          value={segmentedViewValue}
                          onChange={handleSegmentedViewValueChange}
                          data={[
                            { label: 'Table', value: 'grid' },
                            { label: 'Groupings', value: 'group' },
                          ]}
                        />
                        </Grid.Col>
                    </Grid>
                </Grid.Col>
                {segmentedViewValue == "grid" ? gridDimensionOptions : groupHighlightingOptions}
            </Grid>

            <p></p>

            {segmentedViewValue == "grid" ? gridHighlightingOptions : <></>}
        </div>
        

        
        {segmentedViewValue == "grid" ? table : grouping}
        

        </Container>);  // replace with your own

    setVisualization(my_vis_component);

  }, [columnValue, rowValue, altValues, highlightRadioValue, jsonResponses, segmentedViewValue, groupColor]);

  // Grab the LLM(s) response data from the back-end server.
  // Called upon connect to another node, or upon a 'refresh' triggered upstream.
  const grabResponses = () => {
    // For some reason, 'on connect' is called twice upon connection.
    // We detect when an inspector node is already fetching, and disable the second call:
    if (is_fetching) return; 

    // Get the ids from the connected input nodes:
    const input_node_ids = inputEdgesForNode(id).map(e => e.source);

    is_fetching = true;

    // Grab responses associated with those ids:
    fetch(BASE_URL + 'app/grabResponses', {
        method: 'POST',
        headers: {'Content-Type': 'application/json', 'Access-Control-Allow-Origin': '*'},
        body: JSON.stringify({
            'responses': input_node_ids,
        }),
    }).then(function(res) {
        return res.json();
    }).then(function(json) {
        if (json.responses && json.responses.length > 0) {
            setJSONResponses(json.responses);
        }
        is_fetching = false;
    }).catch(() => {
        is_fetching = false; 
    });
  };

  /** Effects to refresh the visualization if anything changes: */
  if (data.input) {
    // If there's a change in inputs to this node, refresh the visualization:
    if (data.input != pastInputs) {
        setPastInputs(data.input);
        grabResponses();
    }
  }

  // Re-grab the responses and recreate the visualization if some other part of the code triggers a 'refresh';
  // for instance, after a prompt node finishes running:
  useEffect(() => {
    if (data.refresh && data.refresh === true) {
        setDataPropsForNode(id, { refresh: false });
        grabResponses();
    }
  }, [data, id, grabResponses, setDataPropsForNode]);

  // The React HTML component to display:
  return (
    <div className="inspector-node cfnode">
    <NodeLabel title={data.title || 'Grid Inspect Node'} 
                nodeId={id}
                icon={'⊞'} />
      <div className="katy-inspect-response-container nowheel nodrag">
        {visualization}
      </div>
      <Handle
        type="target"
        position="left"
        id="input"
        style={{ top: "50%", background: '#555' }}
        onConnect={grabResponses}
      />
    </div>
  );
};

export default GridInspectNode;