import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Handle } from 'react-flow-renderer';
import useStore from './store';
import NodeLabel from './NodeLabelComponent'
import {BASE_URL} from './store';
import { Grid, Select, Radio, NumberInput, Group } from '@mantine/core';
import './output-grid.css';

import embeddingsFilePath from './nlp_models/glove.6B.50d.10k.txt';

import winkNLP from 'wink-nlp'
import winkModel from 'wink-eng-lite-web-model'
import BM25Vectorizer from 'wink-nlp/utilities/bm25-vectorizer'
import similarity from 'wink-nlp/utilities/similarity.js'

const nlp = winkNLP(winkModel);
const its = nlp.its;
const as = nlp.as;


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
  const [sentNum, setSentNum] = useState([]);



  // Update the visualization whenever 'jsonResponses' changes:
  useEffect(() => {
    if (!jsonResponses || (Array.isArray(jsonResponses) && jsonResponses.length === 0))
        return;

    // === Construct a visualization using jsonResponses here ===
    // ....

    

    console.log("similarity of dog and cat", cosine_similarity(embeddings["cat"], embeddings["dog"]));
    console.log("similarity of dog and table", cosine_similarity(embeddings["table"], embeddings["dog"]));

    const badNgram = (ngram) => {
      const cleanNgram = ngram.trim();
      const badWords = [
        "the",
        "and",
        "an",
        "a",
        "as",
        "with",
        "its",
        "of",
        "our",
        "in",
        "this",
        "that",
        "these",
        "those",
        "their",
        "they",
        "them",
        "to",
        ".",
        ",",
        ";"
      ];
      const badBigrams = badWords.flatMap((word1, index1) => 
        badWords.filter((_, index2) => index1 !== index2).map(word2 => word1 + ' ' + word2)
      );
      if (badWords.includes(cleanNgram.toLowerCase())) {
        return true;
      }
      if (badBigrams.includes(cleanNgram.toLowerCase())) {
        return true;
      }
      return false;
    }

    // assumes tokens is a list of list of tokens
    const selectNgrams = (tokens, max_n, k) => {
      let ngrams = {};
      for (let i = 0; i < tokens.length; ++i) {
        for (let j = 0; j < tokens[i].length; ++j) {
          for (let n = 1; n <= max_n; n++) {
            if (j + n > tokens[i].length) {
              continue;
            }
            let ngram = tokens[i].slice(j, j + n);
            ngram = ngram.join(" ");
            if (!ngrams[ngram]) {
              ngrams[ngram] = [];
            }
            ngrams[ngram].push({ i: i, j: j, n: n });
          }
        }
      }

      for (var key in ngrams) {
        if (badNgram(key)) {
          delete ngrams[key];
        }
      }

      // Convert to [key, value] pairs and calculate max n for each ngram
      let pairs = Object.entries(ngrams).map(([ngram, locations]) => {
        let maxN = Math.max(...locations.map((location) => location.n));
        let count = locations.length;
        return [ngram, locations, maxN, count];
      });

      // Sort by max n, in descending order, then by count, also in descending order
      pairs.sort((a, b) => b[2] - a[2] || b[3] - a[3]);

      // Get top k ngrams
      let topKNgrams = {};
      let blockedIndices = Array(tokens.length)
        .fill()
        .map(() => new Set());

      for (let i = 0; i < pairs.length; ++i) {
        const [ngram, locations] = pairs[i];

        if (locations.length === 1) {
          continue;
        }

        // check if an ngram is overlapping with an already-added ngram
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
          topKNgrams[ngram] = locations;
          if (Object.keys(topKNgrams).length === k) break;

          for (const location of locations) {
            for (let idx = location.j; idx < location.j + location.n; idx++) {
              blockedIndices[location.i].add(idx);
            }
          }
        }
      }

      return topKNgrams;
    }

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
       responses is an array of responseTokenized; returns dict of overlapping
       strings and their locations in the format:
       {"text of string": [
            {i: index of response string occurs in,
             j: index of token in response where string starts,
             n: number of tokens in string}
       ]}
    */
    const findOverlap = (responses) => {
        let overlaps = {}
        for (let i = 0; i < responses.length-1; i++) {
            for (let j = i+1; j < responses.length; j++) {
                let lcs = LCSubStr(responses[i], responses[j]);
                if (!overlaps.hasOwnProperty(lcs.string)) {
                    overlaps[lcs.string] = []
                }
                overlaps[lcs.string].push({"i": i, "j": lcs.x, "n": lcs.string.length})
                overlaps[lcs.string].push({"i": j, "j": lcs.y, "n": lcs.string.length})
            }
        }
        return overlaps;
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
            const doc = nlp.readDoc(response.replaceAll("\n", " <br/> "));
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
                bow: nlp.readDoc(sentencesOut[i]).tokens().out(its.value, as.bow)
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

    // learn the bm25 vectors
    const bm25 = BM25Vectorizer();
    const corpus = jsonResponsesMod.map((obj) => obj.response);
    corpus.forEach((doc) => bm25.learn(nlp.readDoc(doc).tokens().out(its.normal)));;

    // pull out tf-idf and top terms
    const idf_array = bm25.out(its.idf);
    const idf_obj = idf_array.reduce((obj, [key, val]) => {
        obj[key] = val;
        return obj;
    }, {});
    console.log("idf_obj", idf_obj);
    const example = bm25.doc(0).out(its.tf);
    const docBowArray = bm25.out(its.docBOWArray);
    // console.log("jsonResponsesMod:", jsonResponsesMod[0].response);
    // console.log("corpus", corpus[0]);
    // console.log("bm25 doc", bm25.doc(0).out(its.normal));
    // console.log("tf", example);
    // console.log("tf-idf", example.map((arr) => [arr[0], arr[1]*idf_obj[arr[0]]]));
    // console.log("docBOWarray", docBowArray);

    const topNTerms = (i, n) => {
        // this version calculate the tf-df for doc number i
        // then returns the top scoring n terms
        const tf = bm25.doc(i).out(its.tf);
        const tfidf = tf.map((arr) => [arr[0], arr[1]*idf_obj[arr[0]]])
        tfidf.sort((a, b) => b[1] - a[1]);
        let top = tfidf.slice(0, n);
        return top.map(term => term[0]);

        // this version takes a docBOWarray (an array of [token, tf] pairs)
        // and returns the top n tokens with the highest tf score
        // return Object.keys(obj).sort((a, b) => obj[b] - obj[a]).slice(0, n)
    };
    const tfidfArray = [...Array(jsonResponsesMod.length).keys()].map((i) => topNTerms(i, 5));
    console.log('tfidfArray', tfidfArray);
    // console.log("docBowArray[0]", docBowArray[0]);
    // console.log("topNTerms", topNTerms(0, 5) );

    // parse sentences and create clusters
    // this is a shitty version. for each sentence in the first response
    // find the n most similar sentences. so there is a "cluster"
    // for each sentence in the first response.
    const coreSentenceSet = jsonResponsesMod[0].sentences;
    const sentenceSimilarities = coreSentenceSet.map((coreSentObj) => {
        // for each "cluster" (i.e. sentence in first response)
        // we will create an array of arrays
        // where each item in the array represents one response
        return jsonResponsesMod.map((respObj) => {
            // for each response we create an array of scores
            // where the scores represent the similarity between
            // that sentence and the core sentence we are looking at
            return respObj.sentences.map((thisSentObj) => {
                return similarity.bow.cosine(coreSentObj.bow, thisSentObj.bow);
            })
        })
    });
    console.log("coreSentenceSet", coreSentenceSet.map((obj) => obj.text));
    console.log("sentenceSimilarities", sentenceSimilarities);
    const bestSentenceSimilarities = sentenceSimilarities.map((cluster) => {
        return cluster.map((scoreList) => {
            // index of max value
            // from https://stackoverflow.com/questions/11301438/return-index-of-greatest-value-in-an-array
            return scoreList.reduce((iMax, x, i, arr) => x > arr[iMax] ? i : iMax, 0);
        })
    });
    console.log("bestSentenceSimilarities", bestSentenceSimilarities);



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
            let cellData = jsonResponsesMod.filter(resp => resp.vars[columnValue] === col_names[x] && resp.vars[rowValue] === row_names[y]);
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
    const rowLCS = gridResponses.map((row) => findOverlap(row.map((cell) => cell.responseTokenized)));
    const colLCS = gridResponses[0].map((_, i) => findOverlap(gridResponses.map((row) => row[i].responseTokenized)));
    const allLCS = findOverlap(gridResponses.flatMap(item => item).map(item => item.responseTokenized));
    console.log("rowLCS", rowLCS);
    console.log("colLCS", colLCS);
    console.log("allLCS", allLCS);


    const rowNgrams = gridResponses.map((row) => selectNgrams(row.map((cell) => cell.responseTokenized), 5, 3));
    const colNgrams = gridResponses[0].map((_, i) => selectNgrams(gridResponses.map((row) => row[i].responseTokenized), 5, 3));
    console.log('rowNgrams', rowNgrams);
    console.log('colNgrams', colNgrams);

    // lcs highlighting uses the same functions as the old ngram ones
    // because the detecting algorithms return the same output form
    const shouldHighlightRowNgrams = (rowIndex, cellIndex, tokenIndex) => {
        for (const [ngram, locations] of Object.entries(rowLCS[rowIndex])) {
            for (let l = 0; l < locations.length; l++) {
              if (locations[l].i === cellIndex) {
                // if this ngram occurs in this string
                if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
                    return true;
                }
              }
            }
        }
        return false;
    }
    const shouldHighlightColNgrams = (rowIndex, cellIndex, tokenIndex) => {
        for (const [ngram, locations] of Object.entries(colLCS[cellIndex])) {
            for (let l = 0; l < locations.length; l++) {
              if (locations[l].i === rowIndex) {
                // if this ngram occurs in this string
                if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
                    return true;
                }
              }
            }
        }
        return false;
    }
    const shouldHighlightAllNgrams = (rowIndex, cellIndex, tokenIndex) => {
        const numCols = gridResponses[0].length;
        for (const [ngram, locations] of Object.entries(allLCS)) {
            for (let l=0; l<locations.length; l++) {
                const thisIndex = (rowIndex * numCols) + cellIndex;
                if (locations[l].i === thisIndex) {
                    if (tokenIndex >= locations[l].j && tokenIndex < locations[l].j + locations[l].n) {
                        return true;
                    }
                }
            }
        }
        return false
    }
    const shouldHighlightTfidf = (cell, tokenIndex) => {
        if (tfidfArray[cell.index].includes(cell.responseTokenized[tokenIndex].toLowerCase())) {
            return true;
        }
        return false;
    }
    const shouldHighlightSentence = (cell, tokenIndex) => {
        if (sentNum == null || sentNum.length == 0) {
            return false;
        }
        const clusterScores = bestSentenceSimilarities[parseInt(sentNum)];
        const bestSentenceIndex = clusterScores[cell.index];
        if (cell.sentences[bestSentenceIndex].tokenIndices.includes(tokenIndex)) {
            return true;
        }
        return false;
    }


    col_names.unshift(""); // add empty string for row headers
    const header = col_names.map((el) => {
            return <th>{el}</th>
        });
    const cells = gridResponses.map((row, rowIndex) => {
        return (
            <tr>
                <th>{row_names[rowIndex]}</th>
                {row.map((cell, cellIndex) => (
                    <td key={cellIndex}>
                            {cell.responseTokenized.map((token, tokenIndex) => {
                                let highlight = false; 
                                let oddEven = false;
                                if (highlightRadioValue === "row") { 
                                    highlight = shouldHighlightRowNgrams(rowIndex, cellIndex, tokenIndex); 
                                    oddEven = rowIndex % 2 == 0;
                                }
                                else if (highlightRadioValue === "col") { 
                                    highlight = shouldHighlightColNgrams(rowIndex, cellIndex, tokenIndex); 
                                    oddEven = cellIndex % 2 == 0;
                                }
                                else if (highlightRadioValue === "lcs") {
                                    highlight = shouldHighlightAllNgrams(rowIndex, cellIndex, tokenIndex);
                                }
                                else if (highlightRadioValue === "tfidf") {
                                    highlight = shouldHighlightTfidf(cell, tokenIndex);
                                }
                                else if (highlightRadioValue === "sent") {
                                    highlight = shouldHighlightSentence(cell, tokenIndex);
                                }
                                const spanStyle = highlight ? {backgroundColor: oddEven ? 'thistle' : 'plum'} : {};
                                if (token == "<br/>") {return <br/>;}
                                if (token == "<br/><br/>") {return <><br/><br/></>;}
                                return <span key={tokenIndex} style={spanStyle}>{token} </span>;;
                            })}
                    </td>
                ))}
            </tr>
        );
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

    const handleSentNumChange = (new_val) => {
        console.log("setting sentNum to", new_val);
        setSentNum(new_val);
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
            <Grid.Col span={2}>
                <Select
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
    const my_vis_component = (<div>
        
        <Grid>
            <Grid.Col span={2}>
            <p>Set table dimensions:</p>
            </Grid.Col>
            <Grid.Col span={2}>
                <Select clearable
                  onChange={handleColumnValueChange}
                  label="Columns"
                  placeholder="Pick one"
                  defaultValue="Model"
                  data={tab_options}
                  value={columnValue}
                />
            </Grid.Col>
            <Grid.Col span={2}>
                <Select clearable
                  onChange={handleRowValueChange}
                  label="Rows"
                  placeholder="Pick one"
                  defaultValue="Model"
                  data={tab_options}
                  value={rowValue}
                />
            </Grid.Col>
            
        </Grid>
        <p></p>
        
        <Grid>
            <Grid.Col span={2}>
                <p>Select which extra variables to view:</p>
            </Grid.Col>
            {alt_select_obj}
        </Grid>
        <p></p>

        <Grid>
            <Grid.Col span={2}>
                <p>Select what to highlight:</p>
            </Grid.Col>
            <Grid.Col span={10}>
                <Radio.Group
                  value={highlightRadioValue}
                  onChange={handleHighlightRadioValue}
                  name="highlightRadioValue"
                  defaultValue="none"
                  // label="Select what you would like to highlight"
                >
                  <Group mt="xs">
                    <Radio value="none" label="None" />
                    <Radio value="row" label="LCS in Rows" />
                    <Radio value="col" label="LCS in Columns" />
                    <Radio value="lcs" label="LCS in all" />
                    <Radio value="tfidf" label="High tf-idf" />
                    <Radio value="sent" label="Similar sentences" />
                  </Group>
                </Radio.Group>
            </Grid.Col>
        </Grid>
        
        <p></p>
        <Grid>
            <Grid.Col span={2}>
                <NumberInput
                  // value="sentNum"
                  onChange={handleSentNumChange}
                  defaultValue={2}
                  min={0}
                  placeholder="Num"
                  label="Sentence number"
                />
            </Grid.Col>
        </Grid>

        <p class="prompt">{prompt}</p>
        
        <table class='outputgrid'>
            <tbody>
                <tr>{header}</tr>
                {cells}
            </tbody>
        </table>
        <div>Hello! The response JSON I received is:
                                <p style={{fontSize: '9pt'}}>{JSON.stringify(jsonResponses)}</p>
                             </div></div>);  // replace with your own
    setVisualization(my_vis_component);

  }, [columnValue, rowValue, altValues, highlightRadioValue, jsonResponses, sentNum]);

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