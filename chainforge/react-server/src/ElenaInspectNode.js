import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Handle } from 'react-flow-renderer';
import useStore from './store';
import NodeLabel from './NodeLabelComponent'
import {BASE_URL} from './store';
import { Grid, Select, Radio, Group } from '@mantine/core';
import './output-grid.css';


const ElenaInspectNode = ({ data, id }) => {

  let is_fetching = false;

  const [visualization, setVisualization] = useState([]);
  const [jsonResponses, setJSONResponses] = useState(null);

  const [pastInputs, setPastInputs] = useState([]);
  const inputEdgesForNode = useStore((state) => state.inputEdgesForNode);
  const setDataPropsForNode = useStore((state) => state.setDataPropsForNode);

  const [columnValue, setColumnValue] = useState([]);
  const [rowValue, setRowValue] = useState([]);
  const [highlightRadioValue, setHighlightRadioValue] = useState([]);

  // Update the visualization whenever 'jsonResponses' changes:
  useEffect(() => {
    if (!jsonResponses || (Array.isArray(jsonResponses) && jsonResponses.length === 0))
        return;

    // === Construct a visualization using jsonResponses here ===
    // ....

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

    const tokenize = (string) => {
        return string.split(" ");
    }

    console.log("CUSTOM INSPECT NODE");
    console.log(jsonResponses);
    // first create an object for each response, copying all other attributes
    // (originally 'response' is a list of responses, but creating one per object makes some other logic simpler)
    // also add 'llm' attribute to 'vars' (because we're going to use 'vars' to set the cols and rows)
    // and create a tokenized version of the response
    let jsonResponsesMod = jsonResponses.flatMap(obj =>
        obj.responses.map((response, index) => ({
            ...obj,
            vars: { ...obj.vars, responseNum: index, llm: obj.llm },
            response: response,
            responseTokenized: tokenize(response)

        }))
    );
    console.log('new', jsonResponsesMod);

    // find all vars for table
    const var_names = [...new Set(jsonResponsesMod.map((obj) => Object.keys(obj.vars)).flat(1))];
    const tab_options = [var_names].flat(1);

    // set the var for rows and volumns
    let col_names = [...new Set(jsonResponsesMod.map((obj) => obj.vars[columnValue]))];
    let row_names = [...new Set(jsonResponsesMod.map((obj) => obj.vars[rowValue]))];

    // look for the data for each cell
    const gridResponses = [];
    for (let y=0; y<row_names.length; y++) {
        let row = [];
        for (let x=0; x<col_names.length; x++) {
            // find the data that goes in that cell
            // assumes it should only find one (i.e. takes first item in filtered list)
            let cellData = jsonResponsesMod.filter(resp => resp.vars[columnValue] === col_names[x] && resp.vars[rowValue] === row_names[y])[0];
            if (!cellData) {
                cellData = {'responseTokenized': ['∅']};
            }

            row.push(cellData);
        }
        gridResponses.push(row);
    }
    console.log('gridResponses', gridResponses);

    const rowNgrams = gridResponses.map((row) => selectNgrams(row.map((cell) => cell.responseTokenized), 5, 3));
    const colNgrams = gridResponses[0].map((_, i) => selectNgrams(gridResponses.map((row) => row[i].responseTokenized), 5, 3));
    console.log('rowNgrams', rowNgrams);
    console.log('colNgrams', colNgrams);

    const shouldHighlightRowNgrams = (rowIndex, cellIndex, tokenIndex) => {
        for (const [ngram, locations] of Object.entries(rowNgrams[rowIndex])) {
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
        for (const [ngram, locations] of Object.entries(colNgrams[cellIndex])) {
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
                                const spanStyle = highlight ? {backgroundColor: oddEven ? 'thistle' : 'plum'} : {};
                                return <span key={tokenIndex} style={spanStyle}>{token} </span>;
                            })}
                    </td>
                ))}
            </tr>
        );
    });

    // ==========================================================

    const handleColumnValueChange = (new_val) => {
        setColumnValue(new_val);
    };
    const handleRowValueChange = (new_val) => {
        setRowValue(new_val);
    };
    const handleHighlightRadioValue = new_val => {
        setHighlightRadioValue(new_val);
    }

    // Set the HTML / React element
    const my_vis_component = (<div>
        <Grid>
            <Grid.Col span={4}>
                <Select
                  onChange={handleColumnValueChange}
                  label="Columns"
                  placeholder="Pick one"
                  defaultValue="Model"
                  data={tab_options}
                  value={columnValue}
                />
            </Grid.Col>
            <Grid.Col span={4}>
                <Select
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
        <Radio.Group
          value={highlightRadioValue}
          onChange={handleHighlightRadioValue}
          name="highlightRadioValue"
          label="Select what you would like to highlight"
        >
          <Group mt="xs">
            <Radio value="none" label="None" />
            <Radio value="row" label="Rows" />
            <Radio value="col" label="Columns" />
          </Group>
        </Radio.Group>
        
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

  }, [columnValue, rowValue, highlightRadioValue, jsonResponses]);

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
    <NodeLabel title={data.title || 'Elena\'s Inspect Node'} 
                nodeId={id}
                icon={'🔍'} />
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

export default ElenaInspectNode;