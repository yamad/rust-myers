// Myers' approximate string matching algorithm
//
// Ref:
// Myers, G. 1999. A Fast Bit-Vector Algorithm for Approximate String
// Matching Based on Dynamic Programming. Journal of ACM 46(3):395--415.
//
// This is the same algorithm that Illumina BaseSpace adapter trimming
// claims to use.
//
//
// P   query sequence of m chars
// T   text of n chars
// k   allowed mismatches
//
// A D.P. matrix (m + 1) x (n + 1) is constructed to solve the problem
//
// A vertical delta column v for the D.P. matrix can be -1, 0, +1 at
// each element, and can be represented by two bit vectors:
//
//    Pv  positive bit-vector, bit set when vertical delta vj is +1
//    Mv  negative bit vector, bit set when vertical delta vj is -1
//        note that when Pv and Mv are both unset, then vertical delta is 0
//
// We also need a bit to represent whether the current characters are matched or not:
//
//    Eq[i,j]
//        1 if p_i == t_j (that is, ith query char matches jth text char), 0 otherwise
//
//
// Consider, each cell in the dynamic programming (D.P.) matrix has:
//
//    inputs: v_in, h_in, and Eq
//   outputs: v_out, h_out
//
// where each v and h is either -1, 0, +1, and Eq is either 0 or 1
// thus, there are 18 possible combinations of the inputs to each cell
//
// Functions Xv and Xh take the inputs to outputs, as follows:
//
//   Xv     = Eq or Mv_in
//   Pv_out = Mh_in or not (Xv or Ph_in)
//   Mv_out = Ph_in and Xv
//
//   Xh     = Eq or Mh_in
//   Ph_out = Mv_in or not (Xh or Pv_in)
//   Mh_out = Pv_in and Xh
//
// where Pv/Mv and Ph/Mh are the bit vectors representing delta_v and delta_h, as above
//
//
// A set of column vectors Peq is precomputed for every possible text
// character s in an alphabet of size sigma, such that:
//
//   Peq[s](i) = (p_i == s)
//
//
// Start:
//   Pv    = 1 (all 1s)
//   Mv    = 0 (all 0s)
//   Score = m (max possible mismatch count)
use std::collections::HashMap;

struct State {
    pv: u64,
    mv: u64,
    score: usize,
}

impl State {
    pub fn new(m: usize) -> Self {
        State {
            pv: (1 << m) - 1,
            mv: 0,
            score: m,
        }
    }
}


fn scan(query: &str, text: &str, k: usize) -> Vec<usize> {
    let mut matches = Vec::new();

    let mut peq = HashMap::with_capacity(6);
    for s in vec!['A', 'T', 'G', 'C', 'N', 'U'] {
        let mut eqs: u64 = 0;
        for (i, p) in query.char_indices() {
            eqs |= ((p == s) as u64) << i;
        }
        peq.insert(s, eqs);
    }

    let m = query.len();


    let mut state = State::new(m);

    for (j, t) in text.char_indices() {
        let o_eq = peq.get(&t);
        match o_eq {
            None => panic!("Alphabet is missing value"),
            Some(eq) => {
                let xv = eq | state.mv;
                let xh = (((eq & state.pv) + state.pv) ^ state.pv) | eq;

                let mut ph = state.mv | !( xh | state.pv );
                let mut mh = state.pv & xh;

                if ph & (1 << (m-1)) > 0 {
                    state.score += 1;
                } else if mh & (1 << (m-1)) > 0 {
                    state.score -= 1;
                }

                ph <<= 1;
                mh <<= 1;
                state.pv = mh | !(xv | ph);
                state.mv = ph & xv;

                if state.score <= k {
                    matches.push(j);
                }
            },
        }
    }

    matches
}


fn main() {
    println!("Hello, world!");
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_one() {
        let text = "ATTCGTGN";
        let query = "TTC";
        assert_eq!(scan(query, text, 0), [3]);
    }
}
