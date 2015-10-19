package main

import "testing"

func Test_map_function(t *testing.T) {
	var tests = []struct {
		s     string
		want  bool
		mymap map[string]string
	}{
		{"hi", true, map[string]string{"hi": ""}},   // "hi" is in map and should therefore be true
		{"bla", false, map[string]string{"hi": ""}}, // "bla" isn't, i.e., false
	}
	for _, c := range tests {
		got := isStrInMap(c.s, c.mymap)
		if got != c.want {
			t.Errorf("Didn't work. put in %q, Wanted %q, got %q", c.s, c.want, got)
		}
	}
}

func Test_getName_function(t *testing.T) {
	var tests = []struct {
		s    []string
		want string
	}{
		{[]string{"previous_element", "dwajdajwklda;Name=hello;don't take me"}, "hello"},
	}
	for _, c := range tests {
		got := getName(c.s)
		if got != c.want {
			t.Errorf("Didn't work, put in %q, wanted %q, got %q", c.s, c.want, got)
		}
	}
}

func Test_countSnps_function(t *testing.T) {
	var tests = []struct {
		f, s []string
		want float64
	}{
		{[]string{"T", "T", "N"}, []string{"T", "T", "N"}, 1.0},
		{[]string{"N", "T", ""}, []string{"T", "T", ""}, 0.5},
		{[]string{"N", "N", ""}, []string{"T", "T", ""}, 0.0},
		{[]string{"", "", ""}, []string{"", "", ""}, 0},
	}
	for _, c := range tests {
		got := compareSNPs(c.f, c.s)
		if got != c.want {
			t.Errorf("Didn't work, put in %q, wanted %v, got %v", c.s, c.want, got)
		}
	}
}

func Test_compareSNPS_LD_function(t *testing.T) {
	var tests = []struct {
		f, s []string
		want float64
	}{
		{[]string{"T", "T", "N"}, []string{"T", "T", "N"}, 1.0},
		{[]string{"T", "T", "T", "N"}, []string{"N", "N", "T", "T"}, 0.3333333333333334},
		{[]string{"", "", ""}, []string{"", "", ""}, 0.0},
	}
	for _, c := range tests {
		got := compareSNPsLD(c.f, c.s)
		if got != c.want {
			t.Errorf("Didn't work, put in %q, wanted %v, got %v", c.s, c.want, got)
		}
	}
}

func Test_compareSNPs_Hamming(t *testing.T) {
	var tests = []struct {
		f, s []string
		want float64
	}{
		{[]string{"T", "T", "N"}, []string{"T", "T", "N"}, 0},
		{[]string{"T", "T", "T", "N"}, []string{"N", "N", "T", "T"}, 3},
		{[]string{"", "", ""}, []string{"", "", ""}, 2.25},
	}
	for _, c := range tests {
		got := compareSNPsHamming(c.f, c.s)
		if got != c.want {
			t.Errorf("Didn't work, put in %q, wanted %v, got %v", c.s, c.want, got)
		}
	}
}
