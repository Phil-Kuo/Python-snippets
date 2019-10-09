def lengthOfLongestSubstring(s):
    """
    :type s: str
    :rtype: int
    """
    if(len(s)==0):
        return 0
    length = len(s)
    max_length = 1
    for i in range(length):
        for j in range(i+1, length+1):
            subset = set(s[i:j])
            tmp = len(subset)
            if((tmp==(j-i)) and (tmp > max_length)):
                max_length = tmp
            else:
                continue
    return max_length


if __name__=="__main__":
    s = "abcabcbb"
    max_length = lengthOfLongestSubstring(s)
    print(max_length)