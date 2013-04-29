<number of states>
<space delimited state names
    (first one is the starting state, last one - ending state;
     there must be no symbol emissions from starting and ending states;
     there must be no transitions to the starting state and from the ending state;
     there must be at least two states: begin and end)
>
<number of different possible symbols (a..z, no more than 26)>
<number of transitions>
<state transitions as space delimited one-per-line triples "from to probability"; unmentioned will have zero probability>
<number of state-symbol emission probablities>
<state-symbol emission probabilities as space delimited one-per-line triples "state symbol probability"; unmentioned will have zero probability>
