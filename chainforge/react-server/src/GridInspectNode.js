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
// https://coolors.co/b4ecb4-f9e2ad-f4c0c0-cfcffc-ade5f9-e3bafd-fcf8b5-fdbaeb
const llmColorPalette = ['#b4ecb4', '#f9e2ad', '#cfcffc', '#ade5f9', '#e3bafd', '#fcf8b5', '#fdbaeb', '#f955cd', '#26e080', '#2654e0', '#7d8191', '#bea5d1'];
// +60 brightness again
// https://coolors.co/e2f8e2-fdf4df-fae5e5-ececfe-dff5fd-f4e3fe-fefce2-fee3f7
const llmColorPaletteDesaturated = ['#e2f8e2', '#fdf4df', '#fae5e5', '#ececfe', '#dff5fd', '#f4e3fe', '#fefce2', '#fee3f7', '#26e080', '#2654e0', '#7d8191', '#bea5d1'];

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
                        if (!sentObj[k].tokenIndices.includes(tokenIndex)) {
                            continue;
                        } else if (sentObj[k].tokenIndices.includes(tokenIndex+numTokens)) {
                            // doesn't cross sentence boundary
                            return [{"i": respIndex, "j": tokenIndex, "n": numTokens, "str": lcsString}];
                        } else {
                            // console.log("crossing", k, tokenIndex, numTokens, lcsString);
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
        });

        // Sort by max n, in descending order, then by count, also in descending order
        pairs.sort((a, b) => b[2] - a[2] || b[3] - a[3]);

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

    // defunct bc of moving to wink-nlp?
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
        corpus.forEach((doc) => bm25.learn(nlp.readDoc(doc).tokens().out(its.normal)));;
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

    const bm25 = getBm25Vectors();
    const idf_obj = getIdfArray(bm25);
    const tfidfArray = [...Array(jsonResponsesMod.length).keys()].map((i) => topTfidfTerms(i, 5));
    console.log('idf_obj', idf_obj);
    console.log('tfidfArray', tfidfArray);


    /*  SIMPLE AGGLOMERATIVE CLUSTERING
        Calculate matrix of similarity scores between all sentences.
        Find most similar two sentences; make them the first cluster.
        Find next most similar two sentences:
            Either add to first cluster, or make them the second cluster.
        Find next most similar two sentences:
            Either add to a cluster, merge two clustesr, or make a new cluster.
        Repeat until n clusters or sentences too similar.
        Return [allSentences, clusters] where:
            allSentences: array of all sentences, where each sentence is
                {bow: {...}, respIndex: int, sentIndex: int, sentence: str }
            clusters: array of clusters, where each cluster is an array of
                sentence ids, i.e. index in allSentences
    */
    const getClusters = () => {
        // First make a list of all sentences
        const allSentences = jsonResponsesMod.map((respObj) => {
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
        // console.log("allSentences", allSentences);

        // then make a similarity matrix
        const n = allSentences.length;
        let simMatrix = Array(n).fill(0).map(() => Array(n).fill(0));
        for (let i=0; i<n; i++) {
            for (let j=0; j<i; j++) {
                const sim = similarity.bow.cosine(allSentences[i].bow, allSentences[j].bow);
                simMatrix[i][j] = sim;
            }
        }
        // console.log("simMatrix", simMatrix);

        // run clustering algo
        let clusters = [];
        const getMaxIndices = (arr) => {
            let max = -Infinity;
            let indices = [];
            for (let i=0; i<arr.length; i++) {
                for (let j=0; j<i; j++) {
                    if (arr[i][j] > max) {
                        max = arr[i][j];
                        indices = [i,j];
                    }
                }
            }
            return [max, indices];
        }
        const checkInCluster = (clusters, i) => {
            for (let k=0; k<clusters.length; k++) {
                if (clusters[k].includes(i)) {
                    return k;
                }
            }
            return null;
        }
        let maxVal; let maxPair;
        [maxVal, maxPair] = getMaxIndices(simMatrix);
        clusters.push(maxPair);
        simMatrix[maxPair[0]][maxPair[1]] = 0;

        const simThreshold = 0.6;
        let clusterPrintOut = []
        while (maxVal > simThreshold) {
            [maxVal, maxPair] = getMaxIndices(simMatrix);

            simMatrix[maxPair[0]][maxPair[1]] = 0;
            // options: add to cluster, merge two clusters, or create new cluster
            let inCluster0 = checkInCluster(clusters, maxPair[0]);
            let inCluster1 = checkInCluster(clusters, maxPair[1]);

            const clusterCopy = clusters.map((cluster) => cluster.slice())
            if (inCluster0 == null && inCluster1 == null) {
                // make new cluster
                clusters.push(maxPair);
            } else if (inCluster0 == null && inCluster1 !== null) {
                // add to cluster
                clusters[inCluster1].push(maxPair[0]);
            } else if (inCluster0 !== null && inCluster1 == null) {
                // add to cluster
                clusters[inCluster0].push(maxPair[1]);
            } else if (inCluster0 !== inCluster1) {
                // merge clusters
                clusters[inCluster0] = clusters[inCluster0].concat(clusters[inCluster1])
                clusters.splice(inCluster1, 1)
            }
        }
        return [allSentences, clusters];
    }

    const [allSentences, clusters] = getClusters();
    console.log("allSentences", allSentences);

    /* SIMPLE SENTENCE GROUPING (EXAMPLORE)
       Calculate the normalized position of each sentence in its response
       Calculate the avg normalized position for each cluster
       Return groupings: ordered array of groups, where each group is an array
        of sentence ids, i.e. index in allSentences
    */
    const getGroupings = () => {

        // function to get normalized position of each sentence
        const calcPosition = (sentId) => {
            const idx = allSentences[sentId].respIndex;
            const numSentences = jsonResponsesMod[idx].sentences.length - 1;
            return allSentences[sentId].sentIndex / numSentences;
        }

        // get the normalized position of each sentence in each cluster
        const positions = clusters.map((cluster) => {
            return cluster.map((sentId) => {
                return calcPosition(sentId);
            });
        });
        console.log("positions", positions);

        // calculate avg position for each cluster
        const avgPosition = positions.map((arr) => {
            const sum = arr.reduce((partialSum, a) => partialSum + a, 0);
            return sum / arr.length;
        });
        console.log("avgPosition", avgPosition);

        // get cluster indices ordered by avg position
        const sortedClusters = avgPosition
            .map((val, index) => ({val, index}))
            .sort((a,b) => a.val - b.val)
            .map(obj => obj.index);
        console.log("sortedClusters", sortedClusters);

        return sortedClusters.map(idx => clusters[idx]);
    }

    for (let k=0; k<clusters.length; k++) {
        const clusterSentences = clusters[k].map((index) => [index, allSentences[index].sentence]);
        console.log("cluster", k, clusterSentences);
    }

    const clusterGroupings = getGroupings();
    console.log("clusterGroupings", clusterGroupings);



    // get the original prompt
    let prompt = jsonResponses[0].prompt
    for (const key in jsonResponses[0].vars) {
        prompt = prompt.replace(jsonResponses[0].vars[key], "{"+key+"}");
    }

    // find all vars for table
    const var_names = [...new Set(jsonResponsesMod.map((obj) => Object.keys(obj.vars)).flat(1))];
    const tab_options = [var_names].flat(1);

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
    const allLCS = findOverlap(jsonResponsesMod.map(item => item.responseTokenized), 5);
    const allLCS_colors = getColors(allLCS);
    // console.log("rowLCS", rowLCS);
    // console.log("colLCS", colLCS);
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
        if (tfidfArray[cell.index].includes(cell.responseTokenized[tokenIndex].toLowerCase())) {
            let colorIndex = cell.index;
            if (colorIndex >= llmColorPalette.length) {
                colorIndex = colorIndex - llmColorPalette.length;
            }
            return [true, llmColorPalette[colorIndex]];
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

        if (clusterId !== null) {
            return [true, llmColorPalette[clusterId]]
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
                                if (highlightRadioValue === "lcs") {
                                    [highlight, color] = shouldHighlightAllNgrams(cell, tokenIndex);
                                }
                                else if (highlightRadioValue === "tfidf") {
                                    [highlight, color] = shouldHighlightTfidf(cell, tokenIndex);
                                }
                                else if (highlightRadioValue === "sent") {
                                    [highlight, color] = shouldHighlightSentence(cell, tokenIndex);
                                }
                                let spanStyle = highlight ? {backgroundColor: oddEven ? 'thistle' : 'plum'} : {};
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

    const grouping = clusterGroupings.map((cluster, clusterIndex) => {
        const clusterSentences = cluster.map((sentId, mappedSentIndex, arr) => {
            const respObj = jsonResponsesMod[allSentences[sentId].respIndex];
            const sentIndex = allSentences[sentId].sentIndex;
            const sentObj = respObj.sentences[sentIndex];
            const tokenIndices = respObj.sentences[sentIndex].tokenIndices;
            const filteredTokenIndices = tokenIndices.filter((tokenIndex) => getToken(sentId, tokenIndex) !== "<br/><br/>");
            
            return (
                <p key={sentId} className="sentenceP">
                {filteredTokenIndices.map((tokenIndex, mappedTokenIndex) => {

                    const token = getToken(sentId, tokenIndex);
                    let spanStyle = {backgroundColor: llmColorPalette[clusterIndex]}

                    if (mappedSentIndex > 0) {
                        const prevSentId = arr[mappedSentIndex-1];
                        const prevSentObj = getSentObj(prevSentId);
                        const prevTokenIndices = prevSentObj.tokenIndices;
                        const filteredPrevTokenIndices = prevTokenIndices.filter((tokenIndex) => getToken(prevSentId, tokenIndex) !== "<br/><br/>");

                        let prevTokenIndex = filteredPrevTokenIndices[mappedTokenIndex];
                        let prevToken = getToken(prevSentId, prevTokenIndex);

                        if (prevToken === token) {
                            spanStyle = {backgroundColor: llmColorPalette[clusterIndex], color: 'grey'}
                        }
                    }
                    if (token == "<br/>") {return <></>;}
                    if (token == "<br/><br/>") {return <></>;}
                    return <span style={spanStyle}>{token} </span>
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

    // Set the HTML / React element
    const my_vis_component = (<Container fluid="true">
        <Grid> 
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
            <Grid.Col span={2} style={{backgroundColor: '#E7F5FF', borderRadius: '5px', margin: '10px'}}>
                <div style={{fontSize: '10pt', color: '#777'}}>Set display type:</div>
                <Grid align="flex-end">
                    <Grid.Col>
                    <label class="mantine-InputWrapper-label mantine-Select-label mantine-jkwmhw" for="mantine-iyl5n76na" id="mantine-iyl5n76na-label" style={{opacity: '0'}}>model</label>
                    <SegmentedControl color="blue"
                      value={segmentedViewValue}
                      onChange={handleSegmentedViewValueChange}
                      data={[
                        { label: 'Grid', value: 'grid' },
                        { label: 'Groupings', value: 'group' },
                      ]}
                    />
                    </Grid.Col>
                </Grid>
            </Grid.Col>
        </Grid>

        <p></p>

        

        <p></p>

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
        
        <p></p>

        <Accordion variant="contained" defaultValue="" chevronPosition="left" chevronSize="15px" style={{margin: '20px'}}>
          <Accordion.Item value="prompt">
            <Accordion.Control style={{fontSize: '10pt', fontWeight: '400 !important'}}>Show prompt:</Accordion.Control>
            <Accordion.Panel style={{fontSize: '10pt', whiteSpace: 'pre-line'}}>{prompt}</Accordion.Panel>
          </Accordion.Item>
        </Accordion>

        
        {segmentedViewValue == "grid" ? table : grouping}
        

        </Container>);  // replace with your own

    setVisualization(my_vis_component);

  }, [columnValue, rowValue, altValues, highlightRadioValue, jsonResponses, segmentedViewValue]);

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