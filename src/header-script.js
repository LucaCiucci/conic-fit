
// @ts-check

// we want to find all the math code enclosed in "$..$" and "$$..$$" and replace it with the
// <i-math> and <tex-block> tags respectively.
// To do this, scan all the text nodes in the document and replace the text nodes with the
// new nodes.
// We also have to exclude the <code> and <pre> tags from this process.

const EXCLUDED_TAGS = ['CODE', 'PRE'];

/**
 * Find all the text nodes in the element and its children
 * 
 * @param {Element} element the element to search
 * @param {function(Text): void} f a function to apply to each text node
 */
function forTextNodesIn(element, f) {
    let children = element.childNodes;
    for (let i = 0; i < children.length; i++) {
        let node = children[i];
        if (node.nodeType === Node.TEXT_NODE) {
            f(/** @type {Text} */(node));
        } else {
            // if instanceof HtmlElement
            if (node.nodeType === Node.ELEMENT_NODE) {
                let e = /** @type {Element} */(node);
                if (EXCLUDED_TAGS.includes(e.tagName)) {
                    continue;
                }
                let childTextNodes = forTextNodesIn(/** @type {Element} */(node), f);
            }
        }
    }
}

function replaceMathWithTags() {
    const textNodes = forTextNodesIn(document.body, (textNode) => {
        let text = textNode.textContent;
        if (!text) {
            return;
        }
        let newText = text.replace(/\$\$(.*?)\$\$/g, '<tex-block>$1</tex-block>');
        newText = newText.replace(/\$(.*?)\$/g, '<i-math>$1</i-math>');
        if (newText !== text) {
            // replace the text node with the new html
            let html = newText;
            let div = document.createElement('div');
            div.innerHTML = html;
            let newNodes = div.childNodes;
            if (!textNode.parentNode) {
                return;
            }
            textNode.parentNode.insertBefore(div, textNode);
            textNode.parentNode.removeChild(textNode);
            for (let j = 0; j < newNodes.length; j++) {
                textNode.parentNode.insertBefore(newNodes[j], div);
            }
            textNode.parentNode.removeChild(div);
        }
    });
}

window.addEventListener('load', replaceMathWithTags);