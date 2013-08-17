import unittest
from biopy import sequtil

class TestSeqUtilFunctions(unittest.TestCase):

    def setUp(self):
        self.seq100 = 'A' * 100

    def test_segment(self):

        # segment sequence of length 100 in 2 segments
        segments = sequtil.segment(self.seq100, 2)

        # check correct functionality
        self.assertEqual(len(segments), 2)
        self.assertTrue(all([len(s) == 50 for s in segments]))

        

if __name__ == '__main__':
    unittest.main()
